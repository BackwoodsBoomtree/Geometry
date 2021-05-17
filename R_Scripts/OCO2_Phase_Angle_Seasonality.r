library(ncdf4)
library(data.table)
library(raster)

file_list   <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2", recursive = TRUE, full.names = TRUE, pattern = "*.nc4")
target_list <- read.csv("C:/Russell/Projects/Geometry/Data/oco2_target_list/OCO2_Target_List.csv")
csv_out     <- "C:/Russell/Projects/Geometry/R_Scripts/CSV/OCO2_PA_Check/OCO2_SIF_vs_PA2.csv"

#### FUNCTIONS ####
cosd <- function(degrees) {
  radians <- cos(degrees * pi / 180)
  return(radians)
}
sind <- function(degrees) {
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
build_data     <- function (input_file) {
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
                   "SIF_757nm" = sif757, "SIF_Uncertainty_757nm" = sif757_U, "continuum_radiance_757nm" = rad757,
                   "lon_1" = lon1, "lon_2" = lon2, "lon_3" = lon3, "lon_4" = lon4,
                   "lat_1" = lat1, "lat_2" = lat2, "lat_3" = lat3, "lat_4" = lat4,
                   "SZA" = sza, "SAz" = saa, "VZA" = vza, "VAz" = vaa, "RAz" = raa, "PA" = pa,
                   "IGBP_index" = igbp)
  df <- na.omit(df) # drop rows that contain an NA anywhere
  # Time Zone
  file_time <<- as.POSIXct((as.numeric(max(df$Delta_Time)) + as.numeric(min(df$Delta_Time))) / 2, origin = "1970-01-01", tz = "UTC")
  file_time <<- as.Date(file_time)
  return(df)
}
subset_flags   <- function(df, mode, flag_cloud, flag_qc) {
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
add_targets    <- function(df, target_list) {
  # Find target
  for (r in 1:nrow(df)) {
    df_distance <- setNames(data.frame(matrix(ncol = 2, nrow = 36)), c("Target_Name", "Distance"))
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
run_regression <- function(df, target_list) {
  
  target_info <- target_list[target_list$Site_Name %like% df$Target_Name[1], ]

  PA  <- abs(df$PA)
  SIF <- df$SIF_757nm
  
  # Model
  model <- lm(SIF ~ PA+I(PA^2))
  myPredict <- predict(model)
  cf <- round(coef(model), 8)
  
  time_diff <- max(df$Delta_Time) - min(df$Delta_Time)
  
  result <- as.data.frame(t(c(target_info$Site_Name, format(df$Delta_Time[1], "%Y-%m-%d"), time_diff, nrow(df), cf[1], cf[2], cf[3],
                              mean(df$SIF_757nm, na.rm = TRUE), target_info$Latitude, target_info$Longitude)))
  
  return(result)
  
}
gen_results    <- function(files, targets, output) {
  
  df_results <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c("Target", "Date", "Time_Diff", "N", "Intercept", "Coeff_1", "Coeff_2", "Mean_SIF757", "Lat", "Lon"))

  for (f in 1:length(files)) {
    df <- build_data(files[f])
    df <- subset_flags(df, 2, 0, 0) # Subset to Target Mode and best quality
    
    if (nrow(df) != 0) {
      print(files[f])
      df <- add_targets(df, targets) # Add target site name to each sounding
      
      for (t in 1:length(unique(df$Target_Name))) {
        group <- df[df$Target_Name == unique(df$Target_Name)[t], ]
        result <- run_regression(group, targets)
        df_results <- rbind(df_results, setNames(result, names(df_results)))
      }
    }
  }
  write.csv(df_results, output, row.names = FALSE)
}


gen_results(file_list, target_list, csv_out)
