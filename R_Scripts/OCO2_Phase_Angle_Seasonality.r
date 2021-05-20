library(ncdf4)
library(data.table)
library(raster)
library(LSD)

file_list      <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2", recursive = TRUE, full.names = TRUE, pattern = "*.nc4")
target_list    <- read.csv("C:/Russell/Projects/Geometry/Data/oco2_target_list/OCO2_Target_List.csv")
csv_out        <- "C:/Russell/Projects/Geometry/R_Scripts/CSV/OCO2_PA_Check/OCO2_SIF_vs_PA_NonABS.csv"
regs_out_dir   <- "C:/Russell/Projects/Geometry/R_Scripts/Figures/OCO2_PA_Check/Regressions"
series_out_dir <- "C:/Russell/Projects/Geometry/R_Scripts/Figures/OCO2_PA_Check/Series"

#### FUNCTIONS ####
cosd                <- function(degrees) {
  radians <- cos(degrees * pi / 180)
  return(radians)
}
sind                <- function(degrees) {
  radians <- sin(degrees * pi / 180)
  return(radians)
}
compute_phase_angle <- function(sat){
  ## Input is satellite data in data.table format with: sza, vza, saa, vaa
  ## Phase angles are set to negative values if the observational azimuth angle is bigger than the solar azimuth angle
  ## (negative phase angle means the sun is to the right of the satellite)
  necessary_names <- c("saa", "sza", "vaa", "vza")
  if (all(necessary_names %in% names(sat))) {
    vaa <- sat$vaa ## viewing azimuth angle
    vza <- sat$vza ## viewing zenith angle
    sza <- sat$sza ## solar zenith angle
    saa <- sat$saa ## solar azimuth angle
    pa  <- phase <- raa <- rep(NA, length(sza))
    phase[vaa > saa] <- -1.
    phase[vaa < saa] <-  1.
    # Relative azimuth angle
    raa <- abs(saa - vaa)
    for (i in 1:length(raa)) {
      if (raa[i] > 180) {
        raa[i] <- abs(raa[i] - 360)
      }
    }
    pa  <- acos(cosd(sza) * cosd(vza) + sind(vza) * sind(sza) * cosd(raa)) * 180. / pi
    pa  <- pa * phase
    return(pa)
  } else {
    print("!!! Necessary input is missing, function returns NULL !!!")
    return(NULL)
  }
}
build_data          <- function (input_file) {
  env <- new.env()
  df <- nc_open(input_file) # Open file
  # Metadata
  mode  <- ncvar_get(df, "Metadata/MeasurementMode") # 0=Nadir, 1=Glint, 2=Target, 3=AreaMap, 4=Transition
  orbit <- ncvar_get(df, "Metadata/OrbitId")
  id    <- ncvar_get(df, "Metadata/SoundingId") # "YYYYMMDDHHMMSS"
  time  <- ncvar_get(df, "Delta_Time") # seconds since 1990-01-01 00:00:00 UTC
  igbp  <- ncvar_get(df, "Science/IGBP_index")
  # SIF and Flag Data
  cloud_flag <- ncvar_get(df, "Cloud/cloud_flag_abp") # 0 - \"Classified clear\", 1 - \"Classified cloudy\", 2 - \"Not classified\", all other values undefined; not used in SIF processing
  q_flag     <- ncvar_get(df, "Quality_Flag") # 0 = best (passes quality control + cloud fraction = 0.0); 1 = good (passes quality control); 2 = bad (failed quality control); -1 = not investigated
  sif757     <- ncvar_get(df, "Science/SIF_757nm")
  sif757_R   <- ncvar_get(df, "Science/SIF_Relative_757nm")
  sif757_U   <- ncvar_get(df, "Science/SIF_Uncertainty_757nm")
  # Radiance
  rad757 <- ncvar_get(df, "Science/continuum_radiance_757nm")
  # Geo Data
  lon_center  <- ncvar_get(df, "Geolocation/longitude")
  lat_center  <- ncvar_get(df, "Geolocation/latitude")
  lon_corners <- ncvar_get(df, "Geolocation/footprint_longitude_vertices")
  lat_corners <- ncvar_get(df, "Geolocation/footprint_latitude_vertices")
  sza <- ncvar_get(df, "SZA")
  saa <- ncvar_get(df, "SAz")
  vza <- ncvar_get(df, "VZA")
  vaa <- ncvar_get(df, "VAz")
  # Close nc file
  nc_close(df)

  # Build Dataframe ####
  # First, get each lat/lon corner into arrays
  lat1 <- lat_corners[1, ]
  lat2 <- lat_corners[2, ]
  lat3 <- lat_corners[3, ]
  lat4 <- lat_corners[4, ]
  lon1 <- lon_corners[1, ]
  lon2 <- lon_corners[2, ]
  lon3 <- lon_corners[3, ]
  lon4 <- lon_corners[4, ]
  # Phase Angle
  pa_table <- data.frame(sza, vza, saa, vaa)  # build table
  pa <- compute_phase_angle(pa_table)
  # Relative azimuth angle
  raa <- abs(saa - vaa)
  for (i in 1:length(raa)) {
    if (raa[i] > 180) {
      raa[i] <- abs(raa[i] - 360)
    }
  }
  df <- data.frame("SoundingID" = id, "MeasurementMode" = mode, "OrbitID" = orbit, "cloud_flag_abp" = cloud_flag, "Quality_Flag" = q_flag,
                   "Delta_Time" = as.POSIXct(time, origin = "1990-01-01", tz = "UTC"), "longitude" = lon_center, "latitude" = lat_center,
                   "SIF_757nm" = sif757, "SIF_Rel_757nm" = sif757_R, "SIF_Uncertainty_757nm" = sif757_U, "continuum_radiance_757nm" = rad757,
                   "lon_1" = lon1, "lon_2" = lon2, "lon_3" = lon3, "lon_4" = lon4,
                   "lat_1" = lat1, "lat_2" = lat2, "lat_3" = lat3, "lat_4" = lat4,
                   "SZA" = sza, "SAz" = saa, "VZA" = vza, "VAz" = vaa, "RAz" = raa, "PA" = pa,
                   "IGBP_index" = igbp)
  df <- na.omit(df) # drop rows that contain an NA anywhere

  return(df)
}
subset_flags        <- function(df, mode, flag_cloud, flag_qc) {
  # mode: 0 = Nadir; 1 = Glint; 2 = Target; 3 = SAM; 4 = Transition; 5 = SAM & Target
  # cloud flag: NA = no filter; 0 = clear; 1 = cloudy; 2 = Not classified; 10 = clear and cloudy
  # qc flag: NA = no filter; 0 = best; 1 = good; 2 = bad; -1 = not investigated; 10 = best and good
  # mode
  if (mode == 0) {
    lab_mode <- "NADIR"
    df <- subset(df, MeasurementMode == 0)
  }else if (mode == 1) {
    lab_mode <- "GLINT"
    df <- subset(df, MeasurementMode == 1)
  } else if (mode == 2) {
    lab_mode <- "TARGET"
    df <- subset(df, MeasurementMode == 2)
  } else if (mode == 3) {
    lab_mode <- "SAM"
    df <- subset(df, MeasurementMode == 3)
  } else if (mode == 4) {
    lab_mode <- "TRANSITION"
    df <- subset(df, MeasurementMode == 4)
  } else if (mode == 5) {
    lab_mode <- "SAM & Target"
    df <- subset(df, MeasurementMode == 2 | MeasurementMode == 3)
  }
  # Cloud
  if (is.na(flag_cloud)) {
    lab_cloud <- "None"
  } else if (flag_cloud == 0) {
    lab_cloud <- "Clear"
    df <- subset(df, cloud_flag_abp == 0)
  } else if (flag_cloud == 1) {
    lab_cloud <- "Cloudy"
    df <- subset(df, cloud_flag_abp == 1)
  } else if (flag_cloud == 2) {
    lab_cloud <- "Not_Classified"
    df <- subset(df, cloud_flag_abp == 2)
  } else if (flag_cloud == 10) {
    lab_cloud <- "Clear?Cloudy"
    df <- subset(df, cloud_flag_abp <= 1)
  }
  # QC
  if (is.na(flag_qc)){
    lab_qc <- "None"
  } else if (flag_qc == 0) {
    lab_qc <- "Best"
    df <- subset(df, Quality_Flag == 0)
  } else if (flag_qc == 1) {
    lab_qc <- "Good"
    df <- subset(df, Quality_Flag == 1)
  } else if (flag_qc == 2) {
    lab_qc <- "Bad"
    df <- subset(df, Quality_Flag == 2)
  } else if (flag_qc == 10) {
    lab_qc <- "Best/Good"
    df <- subset(df, Quality_Flag <= 1)
  } else if (flag_qc == -1) {
    lab_qc <- "Not_Investigated"
    df <- subset(df, Quality_Flag == -1)
  }
  lab_class <- "All" # replaced by subset_cover function if it is used to subset land cover
  list_return <- list("lab_mode" = lab_mode, "lab_cloud" = lab_cloud, "lab_qc" = lab_qc, "lab_class" = lab_class)
  list2env(list_return, .GlobalEnv)
  return(df)
}
add_targets         <- function(df, target_list) {
  # Find target
  for (r in 1:nrow(df)) {
    df_distance <- setNames(data.frame(matrix(ncol = 2, nrow = 28)), c("Target_Name", "Distance"))
    df_distance$Target_Name <- target_list$Site_Name
    
    for (t in 1:nrow(target_list)) {
      # Calculate distance between each sounding and all targets
      # df_distance$Distance[t] <- distm(c(df$longitude[r], df$latitude[r]), c(target_list$Longitude[t], target_list$Latitude[t]), fun = distHaversine)
      df_distance$Distance[t] <- pointDistance(c(df$longitude[r], df$latitude[r]), c(target_list$Longitude[t], target_list$Latitude[t]), lonlat = TRUE)
      
    }
    # Find smallest distance and assign target
    df$Target_Name[r] <- df_distance$Target_Name[df_distance$Distance == min(df_distance$Distance)]
  }
  return(df)
}
plot_regression     <- function(x, y, type, line, prediction, cfs, target_name, date, out_dir) {
  
  par(oma=c(3, 3, 2, 1.25))
  op <- par(mar = c(0, 0, 0, 0))
  cf <- round(cfs, 4) ## rounded coefficients for better output
  ## sign check to avoid having plus followed by minus for negative coefficients
  eq <- paste0("y = ", round(cf[1], 2),
               ifelse(sign(cf[2]) == 1, " + ", " - "), abs(round(cf[2], 2)), "x ",
               ifelse(sign(cf[3]) == 1, " + ", " - "), format(abs(cf[3]), scientific = FALSE), "x^2")
  file_name <- paste0(out_dir, "/", target_name, "_", date, "_", type, ".jpg")
  
  if (type == "SIF") {
    plot_SIF <- function() {
      heatscatter(x, y, cex.axis = 1.25, ylim = c(-1, 2), xlim = c(-80, 80), cexplot = 1, tck = 0.03, mgp = c(3, 0.2, 0),
                  main = "", xlab = NA, ylab = NA)
      mtext(2, text=expression(paste("OCO-2 SIF"['757']*" (mW/m"^{2}*"/sr/nm)")), line = 1.25, cex=1.25)
      if (line == TRUE) { # polynomial line
        ix <- sort(x, index.return = T)$ix
        lines(x[ix], prediction[ix], col = rgb(189, 189, 189, max = 255), lwd = 4)
      }
      title(eq, line = -1) # printing of the equation
      box()
      mtext(3, text = paste0(target_name, " ", date), cex = 1.25)
      mtext(1, text = expression(paste("Phase Angle")), line = 2, cex = 1.25)
    }
    plotit(file_name, sw = 2, sh = 2, sres = 2, plot_SIF, saveit = TRUE, notinR = TRUE, fileformat = "jpeg")
    
  } else if (type == "SIFrel") {
      plot_SIFrel <- function() {
        heatscatter(x, y, cex.axis = 1.25, ylim = c(-0.03, 0.06), xlim = c(-80, 80), cexplot = 1, tck = 0.03, mgp = c(3, 0.2, 0),
                    main = "", xlab = NA, ylab = NA)
        mtext(2, text=expression(paste("OCO-2 Relative SIF"['757']*" (mW/m"^{2}*"/sr/nm)")), line = 1.25, cex=1.25)
        if (line == TRUE) { # polynomial line
          ix <- sort(x, index.return = T)$ix
          lines(x[ix], prediction[ix], col = rgb(189, 189, 189, max = 255), lwd = 4)
        }
        title(eq, line = -1) # printing of the equation
        box()
        mtext(3, text = paste0(target_name, " ", date), cex = 1.25)
        mtext(1, text = expression(paste("Phase Angle")), line = 2, cex = 1.25)
      }
      plotit(file_name, sw = 2, sh = 2, sres = 2, plot_SIFrel, saveit = TRUE, notinR = TRUE, fileformat = "jpeg")
      
  } else if (type == "Rad") {
    plot_Rad <- function() {
      heatscatter(x, y, cex.axis = 1.25, ylim = c(0, 150), xlim = c(-80, 80), cexplot = 1, tck = 0.03, mgp = c(3, 0.2, 0),
                  main = "", xlab = NA, ylab = NA)
      mtext(2, text=expression(paste("OCO-2 Radiance"['757']*" (mW/m"^{2}*"/sr/nm)")), line = 1.25, cex=1.25)
      if (line == TRUE) { # polynomial line
        ix <- sort(x, index.return = T)$ix
        lines(x[ix], prediction[ix], col = rgb(189, 189, 189, max = 255), lwd = 4)
      }
      title(eq, line = -1) # printing of the equation
      box()
      mtext(3, text = paste0(target_name, " ", date), cex=1.25)
      mtext(1, text = expression(paste("Phase Angle")), line=2, cex=1.25)
    }
    plotit(file_name,sw = 2, sh = 2, sres = 2, plot_Rad, saveit = TRUE, notinR = TRUE, fileformat = "jpeg")
  }
}
run_regression      <- function(df, target_list, plot_regs_flag, regs_out_dir, line) {

  PA     <- df$PA
  SIF    <- df$SIF_757nm
  SIFrel <- df$SIF_Rel_757nm
  Rad    <- df$continuum_radiance_757nm
  
  # Regressions
  model_SIF    <- lm(SIF ~ PA+I(PA^2))
  model_SIFrel <- lm(SIFrel ~ PA+I(PA^2))
  model_Rad    <- lm(Rad ~ PA+I(PA^2))
  
  predict_SIF    <- predict(model_SIF)
  predict_SIFrel <- predict(model_SIFrel)
  predict_Rad    <- predict(model_Rad)
  
  cf_SIF    <- round(coef(model_SIF), 8)
  cf_SIFrel <- round(coef(model_SIFrel), 8)
  cf_Rad    <- round(coef(model_Rad), 8)
  
  # Return results and data
  target_info <- target_list[target_list$Site_Name %like% df$Target_Name[1], ]
  target_name <- target_info$Site_Name
  time_diff   <- max(df$Delta_Time) - min(df$Delta_Time) # Length of time it took to take all soundings
  mid_time    <- min(df$Delta_Time) + (time_diff / 2)
  date        <- format(mid_time, "%Y-%m-%d")
  sza_mean    <- mean(df$SZA, na.rm = TRUE)
  
  result <- as.data.frame(t(c(target_name, format(mid_time, c("%Y-%m-%d %H:%M")), time_diff, sza_mean, nrow(df),
                              cf_SIF[1], cf_SIF[2], cf_SIF[3],
                              cf_SIFrel[1], cf_SIFrel[2], cf_SIFrel[3],
                              cf_Rad[1], cf_Rad[2], cf_Rad[3],
                              mean(df$SIF_757nm, na.rm = TRUE),
                              mean(df$SIF_Rel_757nm, na.rm = TRUE),
                              mean(df$continuum_radiance_757nm, na.rm = TRUE),
                              target_info$Latitude, target_info$Longitude)))
  
  # Plot regressions if the flag == TRUE
  if (plot_regs_flag == TRUE && nrow(df) > 1000){
    plot_regression(PA, SIF, "SIF", line, predict_SIF, cf_SIF, target_name, date, regs_out_dir)
    plot_regression(PA, SIFrel, "SIFrel", line, predict_SIFrel, cf_SIFrel, target_name, date, regs_out_dir)
    plot_regression(PA, Rad, "Rad", line, predict_Rad, cf_Rad, target_name, date, regs_out_dir)
  }
  
  return(result)
  
}
gen_results         <- function(files, targets, output, plot_regs_flag, regs_out_dir, line) {
  
  df_results <- setNames(data.frame(matrix(ncol = 19, nrow = 0)),
                         c("Target", "Date", "Time_Diff", "Mean_SZA", "N",
                           "Intercept_SIF", "Coeff_1_SIF", "Coeff_2_SIF",
                           "Intercept_SIFrel", "Coeff_1_SIFrel", "Coeff_2_SIFrel",
                           "Intercept_Rad", "Coeff_1_Rad", "Coeff_2_Rad",
                           "Mean_SIF757", "Mean_SIF757Rel", "Mean_Rad757",
                           "Lat", "Lon"))

  for (f in 1:length(files)) {
    df <- build_data(files[f])
    df <- subset_flags(df, 2, 0, 0) # Subset to Target Mode and best quality
    
    if (nrow(df) != 0) {
      print(files[f])
      df <- add_targets(df, targets) # Add target site name to each sounding
      
      for (t in 1:length(unique(df$Target_Name))) {
        group <- df[df$Target_Name == unique(df$Target_Name)[t], ]
        result <- run_regression(group, targets, plot_regs_flag, regs_out_dir, line)
        df_results <- rbind(df_results, setNames(result, names(df_results)))
      }
    }
  }
  write.csv(df_results, output, row.names = FALSE)
}
plot_results_series <- function(csv, out_dir) {
  
  df     <- read.csv(csv)
  df     <- df[df$N > 1000, ] # Remove days with low number of soundings

  for (t in 1:length(unique(df$Target))) {
    
    target <- unique(df$Target)[t]
    group  <- df[df$Target == target, ]

    ###### SETUP MAIN PLOT ####
    
    pdf(paste0(out_dir, "/", target, ".pdf"), width=7.5, height=4.75, compress=FALSE)
    
    par(oma=c(4,3.75,1.5,3.75))
    
    xlabs <- format(seq(as.Date("2015-01-01"), as.Date("2021-01-01"), by = "quarter"), format = "%Y-%m")
    xlocs <- seq(as.Date("2015-01-01"), as.Date("2021-01-01"), by = "quarter")
    ############## Plot #####################
    
    # Axes margin
    op <- par(mar = c(0,0,0,0))
    
    plot(as.Date(group$Date), group$Coeff_1, type="n", ann = FALSE, axes = FALSE, xlim = as.Date(c("2015-01-01", "2021-01-01")), ylim = c(-0.08, 0.08))
    abline(h = 0)
    rect(as.Date("2015-05-01"), -1.0, as.Date("2015-08-31"), 1.0, col = rgb(229, 229, 229, max = 255), border = NA)
    rect(as.Date("2016-05-01"), -1.0, as.Date("2016-08-31"), 1.0, col = rgb(229, 229, 229, max = 255), border = NA)
    rect(as.Date("2017-05-01"), -1.0, as.Date("2017-08-31"), 1.0, col = rgb(229, 229, 229, max = 255), border = NA)
    rect(as.Date("2018-05-01"), -1.0, as.Date("2018-08-31"), 1.0, col = rgb(229, 229, 229, max = 255), border = NA)
    rect(as.Date("2019-05-01"), -1.0, as.Date("2019-08-31"), 1.0, col = rgb(229, 229, 229, max = 255), border = NA)
    rect(as.Date("2020-05-01"), -1.0, as.Date("2020-08-31"), 1.0, col = rgb(229, 229, 229, max = 255), border = NA)
    abline(h = 0)
    lines(as.Date(group$Date), group$Coeff_1, type="o", col = rgb(67, 147, 195, max = 255), lty = 1, pch = 16, lwd = 2)
    axis(side = 1, labels = xlabs, las = 2, tck = 0.03, mgp = c(3, 0.2, 0), at = xlocs)
    axis(side = 2, las = 1, tck = 0.03, mgp = c(3, 0.3, 0), at = c(seq(from = -0.08, to = 0.08, by = 0.02)))
    
    par(new=T)
    plot(as.Date(group$Date), group$Mean_SIF757, type = "o", col = rgb(174, 1, 126, max = 255), lty = 1, pch = 16, ylim = c(-2, 2),
         xlim =as.Date(c("2015-01-01", "2021-01-01")), ann = FALSE, axes = FALSE, lwd = 2)
    axis(side = 4, las = 1, tck = 0.03, mgp = c(3, 0.3, 0), at = c(seq(from = -2, to = 2, by = 0.5)))
    
    legend(as.Date("2015-01-01"), 2.0, pch = c(16, 16), lty = c(1, 1), horiz = T, lwd = 2,
           bg = "white", col = c(rgb(67, 147, 195, max = 255), rgb(174, 1, 126, max = 255)),
           c("Coefficient", expression(paste("Mean SIF"['757']*""))))
    
    box()
    
    mtext(2 ,text = "Coefficient", line = 2.5, cex = 1.5)
    mtext(3 ,text = target, line = 0.25, cex = 1.5)
    mtext(4, text = expression(paste("Mean SIF"['757']*"")), line = 2.5, cex = 1.5)
    
    dev.off()
  }
}

# Generate the regression results
gen_results(file_list, target_list, csv_out, TRUE, regs_out_dir, FALSE)

# Plot all targets timeseries
plot_results_series(csv_out, series_out_dir)

df <- build_data("C:/Russell/Projects/Geometry/Data/oco2/2016/09/oco2_LtSIF_160913_20210129t133342z.nc4")
df <- subset_flags(df, 2, 0, 0)
df <- add_targets(df, target_list)
