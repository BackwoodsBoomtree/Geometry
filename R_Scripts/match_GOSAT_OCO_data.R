library(raster)
library(rgdal)
library(ncdf4)

options(scipen = 999)

files_gosat <- list.files(path = "C:/Russell/Projects/Geometry/Data/gosat", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_oco2  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)

pair_files <- function(files_gosat, files_oco2) {
  dates_gosat <- vector()
  dates_oco2  <- vector()
    # Get dates from filenames
  for (i in 1:length(files_gosat)){
    date        <- substr(basename(files_gosat[i]), 13, 18)
    dates_gosat <- c(dates_gosat, date)
  }
    for (i in 1:length(files_oco2)){
    date       <- substr(basename(files_oco2[i]), 12, 17)
    dates_oco2 <- c(dates_oco2, date)
  }
    # Combine dates and filenames into dataframes
  files_gosat <- data.frame(dates_gosat, files_gosat)
  files_oco2  <- data.frame(dates_oco2, files_oco2)
    # Pair the filenames
  paired_files <- data.frame(gosat = vector(), oco2 = vector())
  for (i in 1:nrow(files_gosat)){
    for (j in 1:nrow(files_oco2)){
      if (files_gosat$dates_gosat[i] == files_oco2$dates_oco2[j]){
        paired_files[i, 1] <- files_gosat[i, 2]
        paired_files[i, 2] <- files_oco2[j, 2]
      }
    }
  }
    # If condition above is not met, it inserts NA in the paired_files dataframe row
    # Here we remove them and rename the rows, otherwise nrow() returns wrong number
  paired_files <- na.omit(paired_files)
  rownames(paired_files) <- 1:nrow(paired_files)
  return(paired_files)
}
cosd <- function(degrees) {
  radians <- cos(degrees * pi / 180)
  return(radians)
}
sind <- function(degrees) {
  radians <- sin(degrees * pi / 180)
  return(radians)
}
compute_phase_angle <- function(sat ){
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
    raa              <- vaa - saa
    idx <- which(raa < -180.)
    if (length(idx) > 0) raa[idx] <- raa[idx] + 360.
    idx <- which(raa > 180.)
    if (length(idx) > 0) raa[idx] <- raa[idx] - 360.
    raa <- abs(raa)
    pa  <- acos(cosd(sza) * cosd(vza) + sind(vza) * sind(sza) * cosd(raa)) * 180. / pi
    pa  <- pa * phase
    return(pa)
  } else {
    print("!!! Necessary input is missing, function returns NULL !!!")
    return(NULL)
  }
}
build_data <- function (input_file, platform) {
  env <- new.env()
  nc <- nc_open(input_file) # Open file
  # Metadata
  mode  <- ncvar_get(nc, "Metadata/MeasurementMode") # 0=Nadir, 1=Glint, 2=Target, 3=AreaMap, 4=Transition
  orbit <- ncvar_get(nc, "Metadata/OrbitId")
  id    <- ncvar_get(nc, "Metadata/SoundingId") # "YYYYMMDDHHMMSS"
  time  <- ncvar_get(nc, "Delta_Time") # seconds since 1990-01-01 00:00:00 UTC
  # SIF and Flag Data
  cloud_flag <- ncvar_get(nc, "Cloud/cloud_flag_abp") # 0 - \"Classified clear\", 1 - \"Classified cloudy\", 2 - \"Not classified\", all other values undefined; not used in SIF processing
  q_flag     <- ncvar_get(nc, "Quality_Flag") # 0 = best (passes quality control + cloud fraction = 0.0); 1 = good (passes quality control); 2 = bad (failed quality control); -1 = not investigated
  sif740_D   <- ncvar_get(nc, "Daily_SIF_740nm")
  sif757_D   <- ncvar_get(nc, "Daily_SIF_757nm")
  sif771_D   <- ncvar_get(nc, "Daily_SIF_771nm")
  sif740     <- ncvar_get(nc, "SIF_740nm")
  sif757     <- ncvar_get(nc, "Science/SIF_757nm")
  sif771     <- ncvar_get(nc, "Science/SIF_771nm")
  sif740_U   <- ncvar_get(nc, "SIF_Uncertainty_740nm")
  sif757_U   <- ncvar_get(nc, "Science/SIF_Uncertainty_757nm")
  sif771_U   <- ncvar_get(nc, "Science/SIF_Uncertainty_771nm")
  sif757_R   <- ncvar_get(nc, "Science/SIF_Relative_757nm")
  sif771_R   <- ncvar_get(nc, "Science/SIF_Relative_771nm")
  # Radiance
  rad757 <- ncvar_get(nc, "Science/continuum_radiance_757nm")
  rad771 <- ncvar_get(nc, "Science/continuum_radiance_771nm")
  # Geo Data
  lon_center <- ncvar_get(nc, "Geolocation/longitude")
  lat_center <- ncvar_get(nc, "Geolocation/latitude")
  sza        <- ncvar_get(nc, "SZA")
  saa        <- ncvar_get(nc, "SAz")
  vza        <- ncvar_get(nc, "VZA")
  vaa        <- ncvar_get(nc, "VAz")
  # OCO-specific Data
  if (platform == "OCO") {
    igbp        <- ncvar_get(nc, "Science/IGBP_index")
    lon_corners <- ncvar_get(nc, "Geolocation/footprint_longitude_vertices")
    lat_corners <- ncvar_get(nc, "Geolocation/footprint_latitude_vertices")
    temp        <- ncvar_get(nc, "Meteo/temperature_skin")
    # Get each lat/lon corner into arrays
    lat1 <- lat_corners[1, ]
    lat2 <- lat_corners[2, ]
    lat3 <- lat_corners[3, ]
    lat4 <- lat_corners[4, ]
    lon1 <- lon_corners[1, ]
    lon2 <- lon_corners[2, ]
    lon3 <- lon_corners[3, ]
    lon4 <- lon_corners[4, ]
  }
  # Relative azimuth angle
  raa <- abs(saa - vaa)
  for (i in 1:length(raa)) {
    if (raa[i] > 180) {
      raa[i] <- abs(raa[i] - 360)
    }
  }
  # Phase Angle
  pa_table <- data.frame(sza, vza, saa, vaa)  # build table
  pa <- compute_phase_angle(pa_table)

  # Build Dataframe ####
  if (platform == "OCO"){
    df <- data.frame("SoundingID" = id, "MeasurementMode" = mode, "OrbitID" = orbit, "cloud_flag_abp" = cloud_flag, "Quality_Flag" = q_flag,
                    "Delta_Time" = as.POSIXct(time, origin = "1990-01-01", tz = "UTC"), "longitude" = lon_center, "latitude" = lat_center,
                    "SIF_Daily_740nm" = sif740_D, "SIF_Daily_757nm" = sif757_D, "SIF_Daily_771nm" = sif771_D,
                    "SIF_740nm" = sif740, "SIF_Uncertainty_740nm" = sif740_U,
                    "SIF_757nm" = sif757, "SIF_Uncertainty_757nm" = sif757_U, "SIF_Relative_757nm" = sif757_R, "continuum_radiance_757nm" = rad757,
                    "SIF_771nm" = sif771, "SIF_Uncertainty_771nm" = sif771_U, "SIF_Relative_771nm" = sif771_R, "continuum_radiance_771nm" = rad771,
                    "lon_1" = lon1, "lon_2" = lon2, "lon_3" = lon3, "lon_4" = lon4,
                    "lat_1" = lat1, "lat_2" = lat2, "lat_3" = lat3, "lat_4" = lat4,
                    "IGBP_index" = igbp, "Temp" = temp,
                    "SZA" = sza, "SAz" = saa, "VZA" = vza, "VAz" = vaa, "PA" = pa, "RAz" = raa)
  } else {
    df <- data.frame("SoundingID" = id, "MeasurementMode" = mode, "OrbitID" = orbit, "cloud_flag_abp" = cloud_flag,
                "Quality_Flag_P" = q_flag[1, ], "Quality_Flag_S" = q_flag[2, ],
                "Delta_Time" = as.POSIXct(time, origin = "1990-01-01", tz = "UTC"), "longitude" = lon_center, "latitude" = lat_center,
                "SIF_Daily_740nm_P" = sif740_D[1, ], "SIF_Daily_757nm_P" = sif757_D[1, ], "SIF_Daily_771nm_P" = sif771_D[1, ],
                "SIF_Daily_740nm_S" = sif740_D[2, ], "SIF_Daily_757nm_S" = sif757_D[2, ], "SIF_Daily_771nm_S" = sif771_D[2, ],
                "SIF_740nm_P" = sif740[1, ], "SIF_Uncertainty_740nm_P" = sif740_U[1, ],
                "SIF_740nm_S" = sif740[2, ], "SIF_Uncertainty_740nm_S" = sif740_U[2, ],
                "SIF_757nm_P" = sif757[1, ], "SIF_Uncertainty_757nm_P" = sif757_U[1, ], "SIF_Relative_757nm_P" = sif757_R[1, ], "continuum_radiance_757nm_P" = rad757[1, ],
                "SIF_757nm_S" = sif757[2, ], "SIF_Uncertainty_757nm_S" = sif757_U[2, ], "SIF_Relative_757nm_S" = sif757_R[2, ], "continuum_radiance_757nm_S" = rad757[2, ],
                "SIF_771nm_P" = sif771[1, ], "SIF_Uncertainty_771nm_P" = sif771_U[1, ], "SIF_Relative_771nm_P" = sif771_R[1, ], "continuum_radiance_771nm_P" = rad771[1, ],
                "SIF_771nm_S" = sif771[2, ], "SIF_Uncertainty_771nm_S" = sif771_U[2, ], "SIF_Relative_771nm_S" = sif771_R[2, ], "continuum_radiance_771nm_S" = rad771[2, ],
                "SZA" = sza, "SAz" = saa, "VZA" = vza, "VAz" = vaa, "PA" = pa, "RAz" = raa)
  }
  df <- na.omit(df) # drop rows that contain an NA anywhere
  rownames(df) <- 1:nrow(df)

  nc_close(nc)
  return(df)
}
subset_flags <- function(df, mode, flag_cloud, flag_qc, platform) {
  if (platform == "OCO"){
      # Mode OCO
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
      # QC OCO
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
  }
  if (platform == "GOSAT"){
        # Mode GOSAT
      if (mode == 0) {
        lab_mode <- "OB1D"
        df <- subset(df, MeasurementMode == 0)
      } else if (mode == 1) {
        lab_mode <- "OB2D"
        df <- subset(df, MeasurementMode == 1)
      } else if (mode == 1) {
        lab_mode <- "SPOD"
        df <- subset(df, MeasurementMode == 2)
      }
        # QC GOSAT
      if (is.na(flag_qc)){
        lab_qc <- "None"
      } else if (flag_qc == 0) {
        lab_qc <- "Best"
        df <- subset(df, Quality_Flag_P == 0)
        df <- subset(df, Quality_Flag_S == 0)
      } else if (flag_qc == 1) {
        lab_qc <- "Good"
        df <- subset(df, Quality_Flag_P == 1)
        df <- subset(df, Quality_Flag_S == 1)
      } else if (flag_qc == 2) {
        lab_qc <- "Bad"
        df <- subset(df, Quality_Flag_P == 2)
        df <- subset(df, Quality_Flag_S == 2)
      } else if (flag_qc == 10) {
        lab_qc <- "Best/Good"
        df <- subset(df, Quality_Flag_P <= 1)
        df <- subset(df, Quality_Flag_S <= 1)
      } else if (flag_qc == -1) {
        lab_qc <- "Not_Investigated"
        df <- subset(df, Quality_Flag_P == -1)
        df <- subset(df, Quality_Flag_S == -1)
      }
    }

  # Cloud for GOSAT and OCO
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
    lab_cloud <- "Clear_&_Cloudy"
    df <- subset(df, cloud_flag_abp <= 1)
  }

  list_return <- list("lab_mode" = lab_mode, "lab_cloud" = lab_cloud, "lab_qc" = lab_qc)
  list2env(list_return, .GlobalEnv)
  return(df)
}
match_soundings <- function(paired_file_list) {
  matched_gosat_sounding_list <- data.frame()
  matched_oco2_sounding_list  <- data.frame()

  for (i in 1:nrow(paired_file_list)){
    date <- substr(basename(paired_file_list$gosat[i]), 13, 18)
    print(paste0("Working on row number ", i, " and date ", date, "."))

      # Build data into dataframes
    df_gosat <- build_data(paired_file_list$gosat[i], "GOSAT")
    df_oco2  <- build_data(paired_file_list$oco2[i], "OCO")

      # Subset by mode, cloud flag, and qc flag
    df_gosat <- subset_flags(df_gosat, 0, 0, 0, "GOSAT")
    df_oco2  <- subset_flags(df_oco2, 0, 0, 0, "OCO")

    if (nrow(df_gosat) != 0 && nrow(df_oco2) != 0){
        # Transform dfs into SpatialPointsDataFrame
        # For some reason, SpatialPoints function adds lat/lon columns to output, so we deduct it from input data to avoid duplication for OCO2 data,
        # but for gosat we leave it here because the buffer() removes the lon/lat attribute columns
      coords_gosat           <- as.data.frame(cbind(df_gosat$longitude, df_gosat$latitude))
      colnames(coords_gosat) <- c("longitude", "latitude")
      point_df_gosat         <- SpatialPointsDataFrame(coords_gosat, proj4string = CRS("+init=epsg:4326"), data = df_gosat, coords.nrs = c(8, 9))
      coords_oco2            <- as.data.frame(cbind(df_oco2$longitude, df_oco2$latitude))
      colnames(coords_oco2)  <- c("longitude", "latitude")
      point_df_oco2          <- SpatialPointsDataFrame(coords_oco2, proj4string = CRS("+init=epsg:4326"), data = subset(df_oco2, select = -c(longitude, latitude)), coords.nrs = c(7, 8))

        # Do intersections between GOSAT and OCO
      polydf_gosat                <- buffer(point_df_gosat, width = 5000, dissolve = FALSE) # Create SpatialPolygonDataFrame with 5km radius buffer
      matched_gosat_sounding_list <- rbind(matched_gosat_sounding_list, as.data.frame(intersect(polydf_gosat, point_df_oco2)))
      matched_oco2_sounding_list  <- rbind(matched_oco2_sounding_list, as.data.frame(intersect(point_df_oco2, polydf_gosat)))

      print(paste0("In total there are now ", nrow(matched_gosat_sounding_list), " overlapping GOSAT soundings, and ", nrow(matched_oco2_sounding_list), " overlapping OCO-2 soundings."))
    } else {
      print(paste0("Day ", date, " excluded due to empty dataframe. GOSAT nrow = ", nrow(df_gosat), "; OCO2 nrow =  ", nrow(df_oco2), "."))
    }
  }
  return(list(matched_gosat_sounding_list, matched_oco2_sounding_list))
}
remove_non_veg <- function (matched_gosat_sounding_list, matched_oco2_sounding_list) {
    # Removes gosat sounding if a non-veg oco sounding is within the gosat footprint
    # 13 = Urban; 15 = Snow/Ice; 16 = Barren; 17 = Water bodies
  new_matched_gosat_sounding_list <- matched_gosat_sounding_list
  new_matched_oco2_sounding_list  <- matched_oco2_sounding_list
  count <- 0

  for (i in 1:nrow(matched_oco2_sounding_list)) {
    if (matched_oco2_sounding_list$IGBP_index[i] == 13 || matched_oco2_sounding_list$IGBP_index[i] == 15 || matched_oco2_sounding_list$IGBP_index[i] == 16 || matched_oco2_sounding_list$IGBP_index[i] == 17) {
      bad_gosat <- matched_oco2_sounding_list$SoundingID.1[i]

      if (bad_gosat %in% new_matched_oco2_sounding_list$SoundingID.1) { # gosat sounding may have already been removed, so check first
        new_matched_gosat_sounding_list <- subset(new_matched_gosat_sounding_list, SoundingID != bad_gosat)
        new_matched_oco2_sounding_list  <- subset(new_matched_oco2_sounding_list, SoundingID.1 != bad_gosat)
        count <- count + 1
        print(paste0("Removed GOSAT SoundingID: ", bad_gosat))
      }
    }
  }
  print(paste0("Total number of removed GOSAT soundings: ", count))
  return(list(new_matched_gosat_sounding_list, new_matched_oco2_sounding_list))
}
collapse_soundings <- function(matched_gosat_sounding_list, matched_oco2_sounding_list) {
    # Standard error of the mean
  se <- function(x) sd(x) / sqrt(length(x))

  # Subset OCO2 soundings by GOSAT sounding and calculate mean SIF and SE of the mean, then add it to GOSAT sounding list
  for (i in 1:length(unique(matched_oco2_sounding_list$SoundingID.1))) {

    sub_data <- subset(matched_oco2_sounding_list, SoundingID.1 == unique(matched_oco2_sounding_list$SoundingID.1)[i])

    matched_gosat_sounding_list$N_OCO_Soundings[i]             <- nrow(sub_data)
    matched_gosat_sounding_list$Temp[i]                        <- mean(sub_data$Temp)
    # For some reason, doesn't keep as.POSIXct format; have to recast it with default origin of 1970 if to be used directly
    matched_gosat_sounding_list$Mean_OCO_Delta_Time[i]         <- mean(sub_data$Delta_Time)
    matched_gosat_sounding_list$Time_Difference_Secs[i]        <- abs(difftime(matched_gosat_sounding_list$Delta_Time[i], mean(sub_data$Delta_Time), units = "secs"))

    matched_gosat_sounding_list$Mean_OCO_SIF_Daily_740nm[i]    <- mean(sub_data$SIF_Daily_740nm)
    matched_gosat_sounding_list$Mean_OCO_SIF_Daily_757nm[i]    <- mean(sub_data$SIF_Daily_757nm)
    matched_gosat_sounding_list$Mean_OCO_SIF_Daily_771nm[i]    <- mean(sub_data$SIF_Daily_771nm)
    matched_gosat_sounding_list$Mean_OCO_SIF_740nm[i]          <- mean(sub_data$SIF_740nm)
    matched_gosat_sounding_list$Mean_OCO_SIF_757nm[i]          <- mean(sub_data$SIF_757nm)
    matched_gosat_sounding_list$Mean_OCO_SIF_771nm[i]          <- mean(sub_data$SIF_771nm)
    matched_gosat_sounding_list$Mean_OCO_SIF_Relative_757nm[i] <- mean(sub_data$SIF_Relative_757nm)
    matched_gosat_sounding_list$Mean_OCO_SIF_Relative_771nm[i] <- mean(sub_data$SIF_Relative_771nm)
    matched_gosat_sounding_list$Mean_OCO_PA[i]                 <- mean(abs(sub_data$PA))

    matched_gosat_sounding_list$SE_OCO_SIF_Daily_740nm[i]    <- se(sub_data$SIF_Daily_740nm)
    matched_gosat_sounding_list$SE_OCO_SIF_Daily_757nm[i]    <- se(sub_data$SIF_Daily_757nm)
    matched_gosat_sounding_list$SE_OCO_SIF_Daily_771nm[i]    <- se(sub_data$SIF_Daily_771nm)
    matched_gosat_sounding_list$SE_OCO_SIF_740nm[i]          <- se(sub_data$SIF_740nm)
    matched_gosat_sounding_list$SE_OCO_SIF_757nm[i]          <- se(sub_data$SIF_757nm)
    matched_gosat_sounding_list$SE_OCO_SIF_771nm[i]          <- se(sub_data$SIF_771nm)
    matched_gosat_sounding_list$SE_OCO_SIF_Relative_757nm[i] <- se(sub_data$SIF_Relative_757nm)
    matched_gosat_sounding_list$SE_OCO_SIF_Relative_771nm[i] <- se(sub_data$SIF_Relative_771nm)
  }

  return(matched_gosat_sounding_list)
}
build_polyDF <- function(df) {
  #### Create Polygons for each row ####
  for (i in 1:nrow(df)) {
    x <- c(df$lon_1[i], df$lon_2[i], df$lon_3[i], df$lon_4[i], df$lon_1[i])
    y <- c(df$lat_1[i], df$lat_2[i], df$lat_3[i], df$lat_4[i], df$lat_1[i])
    polly <- Polygon(cbind(x, y))
    polly <- Polygons(list(polly), row.names(df[i, ]))
    if (i == 1) {
      pollyLayer <- SpatialPolygons(list(polly), proj4string = CRS("+init=epsg:4326"))
    } else {
      pollyLayer <- SpatialPolygons(c(slot(pollyLayer, "polygons"), list(polly)), proj4string = CRS("+init=epsg:4326"))
    }
  }
  # Assign Data to polygons as a spatial polygon DF data type
  poly_df <- SpatialPolygonsDataFrame(pollyLayer, df)
}


#region ################# RESULTS ##################
paired_file_list <- pair_files(files_gosat, files_oco2)

matched_sounding_list       <- match_soundings(paired_file_list)
matched_gosat_sounding_list <- as.data.frame(matched_sounding_list[1])
matched_oco2_sounding_list  <- as.data.frame(matched_sounding_list[2])

unfiltered_matched_sounding_means <- collapse_soundings(matched_gosat_sounding_list, matched_oco2_sounding_list)

veg_matched_sounding_list   <- remove_non_veg(matched_gosat_sounding_list, matched_oco2_sounding_list)
veg_matched_gosat_sounding_list <- as.data.frame(veg_matched_sounding_list[1])
veg_matched_oco2_sounding_list  <- as.data.frame(veg_matched_sounding_list[2])

veg_matched_sounding_means <- collapse_soundings(veg_matched_gosat_sounding_list, veg_matched_oco2_sounding_list)

time_veg_matched_sounding_means <- subset(veg_matched_sounding_means, Time_Difference_Secs <= 3600)
vza_time_veg_matched_sounding_means <- subset(time_veg_matched_sounding_means, VZA < 5)
n_vza_time_veg_matched_sounding_means <- subset(vza_time_veg_matched_sounding_means, N_OCO_Soundings >= 10)
temp_n_vza_time_veg_matched_sounding_means <- subset(n_vza_time_veg_matched_sounding_means, Temp >= 278.15) # 5 degree C in Kelvin

#endregion

#region ############# PLOTS WITH DIFFERENT FILTERS - GOSAT P #############

### STATS ###

# Mean up the GOSAT polarizations
mean_GOSAT_SIF_740nm_unfiltered          <- unfiltered_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_veg                 <- veg_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_time_veg            <- time_veg_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_vza_time_veg        <- vza_time_veg_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_n_vza_time_veg      <- n_vza_time_veg_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_temp_n_vza_time_veg <- temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_P

# Run linear regressions
reg_SIF_740nm_unfiltered          <- lm(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered)
reg_SIF_740nm_veg                 <- lm(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg)
reg_SIF_740nm_time_veg            <- lm(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg)
reg_SIF_740nm_vza_time_veg        <- lm(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg)
reg_SIF_740nm_n_vza_time_veg      <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg)
reg_SIF_740nm_temp_n_vza_time_veg <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg)


# Summaries
summary_unfiltered          <- summary(reg_SIF_740nm_unfiltered)
summary_veg                 <- summary(reg_SIF_740nm_veg)
summary_time_veg            <- summary(reg_SIF_740nm_time_veg)
summary_vza_time_veg        <- summary(reg_SIF_740nm_vza_time_veg)
summary_n_vza_time_veg      <- summary(reg_SIF_740nm_n_vza_time_veg)
summary_temp_n_vza_time_veg <- summary(reg_SIF_740nm_temp_n_vza_time_veg)

# Round coefficients
cf_reg_SIF_740nm_unfiltered           <- round(coef(reg_SIF_740nm_unfiltered), 2)
cf_reg_SIF_740nm_veg                  <- round(coef(reg_SIF_740nm_veg), 2)
cf_reg_SIF_740nm_time_veg             <- round(coef(reg_SIF_740nm_time_veg), 2)
cf_reg_SIF_740nm_vza_time_veg         <- round(coef(reg_SIF_740nm_vza_time_veg), 2)
cf_reg_SIF_740nm_n_vza_time_veg       <- round(coef(reg_SIF_740nm_n_vza_time_veg), 2)
cf_reg_SIF_740nm_temp_n_vza_time_veg  <- round(coef(reg_SIF_740nm_temp_n_vza_time_veg), 2)

# Equations
eq_unfiltered     <- paste0("y = ", cf_reg_SIF_740nm_unfiltered[1],
             ifelse(sign(cf_reg_SIF_740nm_unfiltered[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_unfiltered[2]), "x")
eq_veg            <- paste0("y = ", cf_reg_SIF_740nm_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_veg[2]), "x")
eq_time_veg       <- paste0("y = ", cf_reg_SIF_740nm_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_time_veg[2]), "x")
eq_vza_time_veg   <- paste0("y = ", cf_reg_SIF_740nm_vza_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_vza_time_veg[2]), "x")
eq_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_n_vza_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_n_vza_time_veg[2]), "x")
eq_temp_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_temp_n_vza_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]), "x")

# RMSE
RMSE_unfiltered          <- paste0("RMSE = ", round(sqrt(mean(summary_unfiltered$residuals ^ 2)), 2))
RMSE_veg                 <- paste0("RMSE = ", round(sqrt(mean(summary_veg$residuals ^ 2)), 2))
RMSE_time_veg            <- paste0("RMSE = ", round(sqrt(mean(summary_time_veg$residuals ^ 2)), 2))
RMSE_vza_time_veg        <- paste0("RMSE = ", round(sqrt(mean(summary_vza_time_veg$residuals ^ 2)), 2))
RMSE_n_vza_time_veg      <- paste0("RMSE = ", round(sqrt(mean(summary_n_vza_time_veg$residuals ^ 2)), 2))
RMSE_temp_n_vza_time_veg <- paste0("RMSE = ", round(sqrt(mean(summary_temp_n_vza_time_veg$residuals ^ 2)), 2))

# R values
r_unfiltered          <- paste0("p < .001, R2 = ", round(summary_unfiltered$adj.r.squared, 2))
r_veg                 <- paste0("p < .001, R2 = ", round(summary_veg$adj.r.squared, 2))
r_time_veg            <- paste0("p < .001, R2 = ", round(summary_time_veg$adj.r.squared, 2))
r_vza_time_veg        <- paste0("p < .001, R2 = ", round(summary_vza_time_veg$adj.r.squared, 2))
r_n_vza_time_veg      <- paste0("p < .001, R2 = ", round(summary_n_vza_time_veg$adj.r.squared, 2))
r_temp_n_vza_time_veg <- paste0("p < .001, R2 = ", round(summary_temp_n_vza_time_veg$adj.r.squared, 2))

##### SETUP MAIN PLOT #####

# pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_GOSAT_OCO2_again.pdf", width=8, height=8, compress=FALSE)

par(mfrow=c(2, 3), oma = c(0, 0.5, 0, 1.5))

# FILTERED DATA PLOTS

# Unfiltered
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_unfiltered, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_unfiltered, r_unfiltered, RMSE_unfiltered), bty = 'n')
mtext(expression(paste("Unfiltered")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT P SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_veg, r_veg, RMSE_veg), bty = 'n')
mtext(expression(paste("Vegetation Only")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT P SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_time_veg, r_time_veg, RMSE_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT P SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_vza_time_veg, r_vza_time_veg, RMSE_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT P SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_n_vza_time_veg, r_n_vza_time_veg, RMSE_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT P SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N Temp
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_temp_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_temp_n_vza_time_veg, r_temp_n_vza_time_veg, RMSE_temp_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10 / T >= 5C")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT P SIF 740nm")), 1, 1.5, cex = 1)
box()

# dev.off()

#endregion

#region ############# PLOTS WITH DIFFERENT FILTERS - GOSAT S #############

### STATS ###

# Mean up the GOSAT polarizations
mean_GOSAT_SIF_740nm_unfiltered          <- unfiltered_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_veg                 <- veg_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_time_veg            <- time_veg_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_vza_time_veg        <- vza_time_veg_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_n_vza_time_veg      <- n_vza_time_veg_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_temp_n_vza_time_veg <- temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_S

# Run linear regressions
reg_SIF_740nm_unfiltered          <- lm(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered)
reg_SIF_740nm_veg                 <- lm(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg)
reg_SIF_740nm_time_veg            <- lm(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg)
reg_SIF_740nm_vza_time_veg        <- lm(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg)
reg_SIF_740nm_n_vza_time_veg      <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg)
reg_SIF_740nm_temp_n_vza_time_veg <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg)


# Summaries
summary_unfiltered          <- summary(reg_SIF_740nm_unfiltered)
summary_veg                 <- summary(reg_SIF_740nm_veg)
summary_time_veg            <- summary(reg_SIF_740nm_time_veg)
summary_vza_time_veg        <- summary(reg_SIF_740nm_vza_time_veg)
summary_n_vza_time_veg      <- summary(reg_SIF_740nm_n_vza_time_veg)
summary_temp_n_vza_time_veg <- summary(reg_SIF_740nm_temp_n_vza_time_veg)

# Round coefficients
cf_reg_SIF_740nm_unfiltered           <- round(coef(reg_SIF_740nm_unfiltered), 2)
cf_reg_SIF_740nm_veg                  <- round(coef(reg_SIF_740nm_veg), 2)
cf_reg_SIF_740nm_time_veg             <- round(coef(reg_SIF_740nm_time_veg), 2)
cf_reg_SIF_740nm_vza_time_veg         <- round(coef(reg_SIF_740nm_vza_time_veg), 2)
cf_reg_SIF_740nm_n_vza_time_veg       <- round(coef(reg_SIF_740nm_n_vza_time_veg), 2)
cf_reg_SIF_740nm_temp_n_vza_time_veg  <- round(coef(reg_SIF_740nm_temp_n_vza_time_veg), 2)

# Equations
eq_unfiltered     <- paste0("y = ", cf_reg_SIF_740nm_unfiltered[1],
             ifelse(sign(cf_reg_SIF_740nm_unfiltered[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_unfiltered[2]), "x")
eq_veg            <- paste0("y = ", cf_reg_SIF_740nm_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_veg[2]), "x")
eq_time_veg       <- paste0("y = ", cf_reg_SIF_740nm_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_time_veg[2]), "x")
eq_vza_time_veg   <- paste0("y = ", cf_reg_SIF_740nm_vza_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_vza_time_veg[2]), "x")
eq_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_n_vza_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_n_vza_time_veg[2]), "x")
eq_temp_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_temp_n_vza_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]), "x")

# RMSE
RMSE_unfiltered          <- paste0("RMSE = ", round(sqrt(mean(summary_unfiltered$residuals ^ 2)), 2))
RMSE_veg                 <- paste0("RMSE = ", round(sqrt(mean(summary_veg$residuals ^ 2)), 2))
RMSE_time_veg            <- paste0("RMSE = ", round(sqrt(mean(summary_time_veg$residuals ^ 2)), 2))
RMSE_vza_time_veg        <- paste0("RMSE = ", round(sqrt(mean(summary_vza_time_veg$residuals ^ 2)), 2))
RMSE_n_vza_time_veg      <- paste0("RMSE = ", round(sqrt(mean(summary_n_vza_time_veg$residuals ^ 2)), 2))
RMSE_temp_n_vza_time_veg <- paste0("RMSE = ", round(sqrt(mean(summary_temp_n_vza_time_veg$residuals ^ 2)), 2))

# R values
r_unfiltered          <- paste0("p < .001, R2 = ", round(summary_unfiltered$adj.r.squared, 2))
r_veg                 <- paste0("p < .001, R2 = ", round(summary_veg$adj.r.squared, 2))
r_time_veg            <- paste0("p < .001, R2 = ", round(summary_time_veg$adj.r.squared, 2))
r_vza_time_veg        <- paste0("p < .001, R2 = ", round(summary_vza_time_veg$adj.r.squared, 2))
r_n_vza_time_veg      <- paste0("p < .001, R2 = ", round(summary_n_vza_time_veg$adj.r.squared, 2))
r_temp_n_vza_time_veg <- paste0("p < .001, R2 = ", round(summary_temp_n_vza_time_veg$adj.r.squared, 2))

##### SETUP MAIN PLOT #####

# pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_GOSAT_OCO2_again.pdf", width=8, height=8, compress=FALSE)

par(mfrow=c(2, 3), oma = c(0, 0.5, 0, 1.5))

# FILTERED DATA PLOTS

# Unfiltered
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_unfiltered, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_unfiltered, r_unfiltered, RMSE_unfiltered), bty = 'n')
mtext(expression(paste("Unfiltered")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_veg, r_veg, RMSE_veg), bty = 'n')
mtext(expression(paste("Vegetation Only")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_time_veg, r_time_veg, RMSE_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_vza_time_veg, r_vza_time_veg, RMSE_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_n_vza_time_veg, r_n_vza_time_veg, RMSE_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N Temp
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_temp_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_temp_n_vza_time_veg, r_temp_n_vza_time_veg, RMSE_temp_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10 / T >= 5C")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# dev.off()

#endregion

#region ############# PLOTS WITH DIFFERENT FILTERS - GOSAT Mean PS #############

### STATS ###

# Mean up the GOSAT polarizations
mean_GOSAT_SIF_740nm_unfiltered          <- (unfiltered_matched_sounding_means$SIF_740nm_P + unfiltered_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_veg                 <- (veg_matched_sounding_means$SIF_740nm_P + veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_time_veg            <- (time_veg_matched_sounding_means$SIF_740nm_P + time_veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_vza_time_veg        <- (vza_time_veg_matched_sounding_means$SIF_740nm_P + vza_time_veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_n_vza_time_veg      <- (n_vza_time_veg_matched_sounding_means$SIF_740nm_P + n_vza_time_veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_temp_n_vza_time_veg <- (temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_S) / 2

# Run linear regressions
reg_SIF_740nm_unfiltered          <- lm(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered)
reg_SIF_740nm_veg                 <- lm(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg)
reg_SIF_740nm_time_veg            <- lm(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg)
reg_SIF_740nm_vza_time_veg        <- lm(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg)
reg_SIF_740nm_n_vza_time_veg      <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg)
reg_SIF_740nm_temp_n_vza_time_veg <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg)


# Summaries
summary_unfiltered          <- summary(reg_SIF_740nm_unfiltered)
summary_veg                 <- summary(reg_SIF_740nm_veg)
summary_time_veg            <- summary(reg_SIF_740nm_time_veg)
summary_vza_time_veg        <- summary(reg_SIF_740nm_vza_time_veg)
summary_n_vza_time_veg      <- summary(reg_SIF_740nm_n_vza_time_veg)
summary_temp_n_vza_time_veg <- summary(reg_SIF_740nm_temp_n_vza_time_veg)

# Round coefficients
cf_reg_SIF_740nm_unfiltered           <- round(coef(reg_SIF_740nm_unfiltered), 2)
cf_reg_SIF_740nm_veg                  <- round(coef(reg_SIF_740nm_veg), 2)
cf_reg_SIF_740nm_time_veg             <- round(coef(reg_SIF_740nm_time_veg), 2)
cf_reg_SIF_740nm_vza_time_veg         <- round(coef(reg_SIF_740nm_vza_time_veg), 2)
cf_reg_SIF_740nm_n_vza_time_veg       <- round(coef(reg_SIF_740nm_n_vza_time_veg), 2)
cf_reg_SIF_740nm_temp_n_vza_time_veg  <- round(coef(reg_SIF_740nm_temp_n_vza_time_veg), 2)

# Equations
eq_unfiltered     <- paste0("y = ", cf_reg_SIF_740nm_unfiltered[1],
             ifelse(sign(cf_reg_SIF_740nm_unfiltered[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_unfiltered[2]), "x")
eq_veg            <- paste0("y = ", cf_reg_SIF_740nm_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_veg[2]), "x")
eq_time_veg       <- paste0("y = ", cf_reg_SIF_740nm_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_time_veg[2]), "x")
eq_vza_time_veg   <- paste0("y = ", cf_reg_SIF_740nm_vza_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_vza_time_veg[2]), "x")
eq_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_n_vza_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_n_vza_time_veg[2]), "x")
eq_temp_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_temp_n_vza_time_veg[1],
             ifelse(sign(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]), "x")

# RMSE
RMSE_unfiltered          <- paste0("RMSE = ", round(sqrt(mean(summary_unfiltered$residuals ^ 2)), 2))
RMSE_veg                 <- paste0("RMSE = ", round(sqrt(mean(summary_veg$residuals ^ 2)), 2))
RMSE_time_veg            <- paste0("RMSE = ", round(sqrt(mean(summary_time_veg$residuals ^ 2)), 2))
RMSE_vza_time_veg        <- paste0("RMSE = ", round(sqrt(mean(summary_vza_time_veg$residuals ^ 2)), 2))
RMSE_n_vza_time_veg      <- paste0("RMSE = ", round(sqrt(mean(summary_n_vza_time_veg$residuals ^ 2)), 2))
RMSE_temp_n_vza_time_veg <- paste0("RMSE = ", round(sqrt(mean(summary_temp_n_vza_time_veg$residuals ^ 2)), 2))

# R values
r_unfiltered          <- paste0("p < .001, R2 = ", round(summary_unfiltered$adj.r.squared, 2))
r_veg                 <- paste0("p < .001, R2 = ", round(summary_veg$adj.r.squared, 2))
r_time_veg            <- paste0("p < .001, R2 = ", round(summary_time_veg$adj.r.squared, 2))
r_vza_time_veg        <- paste0("p < .001, R2 = ", round(summary_vza_time_veg$adj.r.squared, 2))
r_n_vza_time_veg      <- paste0("p < .001, R2 = ", round(summary_n_vza_time_veg$adj.r.squared, 2))
r_temp_n_vza_time_veg <- paste0("p < .001, R2 = ", round(summary_temp_n_vza_time_veg$adj.r.squared, 2))

##### SETUP MAIN PLOT #####

# pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_GOSAT_OCO2_again.pdf", width=8, height=8, compress=FALSE)

par(mfrow=c(2, 3), oma = c(0, 0.5, 0, 1.5))

# FILTERED DATA PLOTS

# Unfiltered
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_unfiltered, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_unfiltered, r_unfiltered, RMSE_unfiltered), bty = 'n')
mtext(expression(paste("Unfiltered")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_veg, r_veg, RMSE_veg), bty = 'n')
mtext(expression(paste("Vegetation Only")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_time_veg, r_time_veg, RMSE_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_vza_time_veg, r_vza_time_veg, RMSE_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_n_vza_time_veg, r_n_vza_time_veg, RMSE_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N Temp
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_temp_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_temp_n_vza_time_veg, r_temp_n_vza_time_veg, RMSE_temp_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10 / T >= 5C")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# dev.off()

#endregion

#region ############# SIF, DAILY, RELATIVE ################

### STATS ###

# Mean up the GOSAT polarizations
mean_GOSAT_SIF_740nm          <- (n_vza_time_veg_matched_sounding_means$SIF_740nm_P + n_vza_time_veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_771nm          <- (n_vza_time_veg_matched_sounding_means$SIF_757nm_P + n_vza_time_veg_matched_sounding_means$SIF_757nm_S) / 2
mean_GOSAT_SIF_771nm          <- (n_vza_time_veg_matched_sounding_means$SIF_771nm_P + n_vza_time_veg_matched_sounding_means$SIF_771nm_S) / 2
mean_GOSAT_SIF_Relative_757nm <- (n_vza_time_veg_matched_sounding_means$SIF_Relative_757nm_P + n_vza_time_veg_matched_sounding_means$SIF_Relative_757nm_S) / 2
mean_GOSAT_SIF_Relative_771nm <- (n_vza_time_veg_matched_sounding_means$SIF_Relative_771nm_P + n_vza_time_veg_matched_sounding_means$SIF_Relative_771nm_S) / 2
mean_GOSAT_SIF_Daily_740nm    <- (n_vza_time_veg_matched_sounding_means$SIF_Daily_740nm_P + n_vza_time_veg_matched_sounding_means$SIF_Daily_740nm_S) / 2
mean_GOSAT_SIF_Daily_757nm    <- (n_vza_time_veg_matched_sounding_means$SIF_Daily_757nm_P + n_vza_time_veg_matched_sounding_means$SIF_Daily_757nm_S) / 2
mean_GOSAT_SIF_Daily_771nm    <- (n_vza_time_veg_matched_sounding_means$SIF_Daily_771nm_P + n_vza_time_veg_matched_sounding_means$SIF_Daily_771nm_S) / 2

# Run linear regressions
reg_SIF_740nm          <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm)
reg_SIF_757nm          <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_757nm ~ mean_GOSAT_SIF_757nm)
reg_SIF_771nm          <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_771nm ~ mean_GOSAT_SIF_771nm)
reg_SIF_Relative_757nm <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_757nm ~ mean_GOSAT_SIF_Relative_757nm)
reg_SIF_Relative_771nm <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_771nm ~ mean_GOSAT_SIF_Relative_771nm)
reg_SIF_Daily_740nm    <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_740nm ~ mean_GOSAT_SIF_Daily_740nm)
reg_SIF_Daily_757nm    <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_757nm ~ mean_GOSAT_SIF_Daily_757nm)
reg_SIF_Daily_771nm    <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_771nm ~ mean_GOSAT_SIF_Daily_771nm)

# Summaries
summary_SIF_740nm          <- summary(reg_SIF_740nm)
summary_SIF_757nm          <- summary(reg_SIF_757nm)
summary_SIF_771nm          <- summary(reg_SIF_771nm)
summary_SIF_Relative_757nm <- summary(reg_SIF_Relative_757nm)
summary_SIF_Relative_771nm <- summary(reg_SIF_Relative_771nm)
summary_SIF_Daily_740nm    <- summary(reg_SIF_Daily_740nm)
summary_SIF_Daily_757nm    <- summary(reg_SIF_Daily_757nm)
summary_SIF_Daily_771nm    <- summary(reg_SIF_Daily_771nm)

# Round coefficients
cf_reg_SIF_740nm <- round(coef(reg_SIF_740nm), 2)
cf_reg_SIF_757nm <- round(coef(reg_SIF_757nm), 2)
cf_reg_SIF_771nm <- round(coef(reg_SIF_771nm), 2)
cf_reg_SIF_Relative_757nm <- round(coef(reg_SIF_Relative_757nm), 2)
cf_reg_SIF_Relative_771nm <- round(coef(reg_SIF_Relative_771nm), 2)
cf_reg_SIF_Daily_740nm <- round(coef(reg_SIF_Daily_740nm), 2)
cf_reg_SIF_Daily_757nm <- round(coef(reg_SIF_Daily_757nm), 2)
cf_reg_SIF_Daily_771nm <- round(coef(reg_SIF_Daily_771nm), 2)

# Equations
eq_SIF_740nm          <- paste0("y = ", cf_reg_SIF_740nm[1],
             ifelse(sign(cf_reg_SIF_740nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm[2]), "x")
eq_SIF_757nm          <- paste0("y = ", cf_reg_SIF_757nm[1],
             ifelse(sign(cf_reg_SIF_757nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_757nm[2]), "x")
eq_SIF_771nm          <- paste0("y = ", cf_reg_SIF_771nm[1],
             ifelse(sign(cf_reg_SIF_771nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_771nm[2]), "x")
eq_SIF_Relative_757nm <- paste0("y = ", cf_reg_SIF_Relative_757nm[1],
             ifelse(sign(cf_reg_SIF_757nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Relative_757nm[2]), "x")
eq_SIF_Relative_771nm <- paste0("y = ", cf_reg_SIF_Relative_771nm[1],
             ifelse(sign(cf_reg_SIF_771nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Relative_771nm[2]), "x")
eq_SIF_Daily_740nm    <- paste0("y = ", cf_reg_SIF_Daily_740nm[1],
             ifelse(sign(cf_reg_SIF_Daily_740nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Daily_740nm[2]), "x")
eq_SIF_Daily_757nm    <- paste0("y = ", cf_reg_SIF_Daily_757nm[1],
             ifelse(sign(cf_reg_SIF_Daily_757nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Daily_757nm[2]), "x")
eq_SIF_Daily_771nm    <- paste0("y = ", cf_reg_SIF_Daily_771nm[1],
             ifelse(sign(cf_reg_SIF_Daily_771nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Daily_771nm[2]), "x")

# RMSE
RMSE_SIF_740nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_740nm$residuals ^ 2)), 2))
RMSE_SIF_757nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_757nm$residuals ^ 2)), 2))
RMSE_SIF_771nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_771nm$residuals ^ 2)), 2))
RMSE_SIF_Relative_757nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Relative_757nm$residuals ^ 2)), 2))
RMSE_SIF_Relative_771nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Relative_771nm$residuals ^ 2)), 2))
RMSE_SIF_Daily_740nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Daily_740nm$residuals ^ 2)), 2))
RMSE_SIF_Daily_757nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Daily_757nm$residuals ^ 2)), 2))
RMSE_SIF_Daily_771nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Daily_771nm$residuals ^ 2)), 2))

# R values
r_SIF_740nm <- paste0("p < .001, R2 = ", round(summary_SIF_740nm$adj.r.squared, 2))
r_SIF_757nm <- paste0("p < .001, R2 = ", round(summary_SIF_757nm$adj.r.squared, 2))
r_SIF_771nm <- paste0("p < .001, R2 = ", round(summary_SIF_771nm$adj.r.squared, 2))
r_SIF_Relative_757nm <- paste0("p < .001, R2 = ", round(summary_SIF_Relative_757nm$adj.r.squared, 2))
r_SIF_Relative_771nm <- paste0("p < .001, R2 = ", round(summary_SIF_Relative_771nm$adj.r.squared, 2))
r_SIF_Daily_740nm <- paste0("p < .001, R2 = ", round(summary_SIF_Daily_740nm$adj.r.squared, 2))
r_SIF_Daily_757nm <- paste0("p < .001, R2 = ", round(summary_SIF_Daily_757nm$adj.r.squared, 2))
r_SIF_Daily_771nm <- paste0("p < .001, R2 = ", round(summary_SIF_Daily_771nm$adj.r.squared, 2))

##### SETUP MAIN PLOT #####

# pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_GOSAT_OCO2_again.pdf", width=8, height=8, compress=FALSE)

par(mfrow=c(3, 3), oma = c(0, 0.5, 0.5, 0.5))

# FILTERED DATA PLOTS

# SIF 757
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_757nm ~ mean_GOSAT_SIF_757nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_757nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_757nm, r_SIF_757nm, RMSE_SIF_757nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF 757nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 757nm")), 1, 1.5, cex = 1)
box()

# SIF 771
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_771nm ~ mean_GOSAT_SIF_771nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_771nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_771nm, r_SIF_771nm, RMSE_SIF_771nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF 771nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 771nm")), 1, 1.5, cex = 1)
box()

# SIF 740
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_740nm, r_SIF_740nm, RMSE_SIF_740nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# SIF DAILY 757
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_757nm ~ mean_GOSAT_SIF_Daily_757nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Daily_757nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Daily_757nm, r_SIF_Daily_757nm, RMSE_SIF_Daily_757nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Daily 757nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Daily 757nm")), 1, 1.5, cex = 1)
box()

# SIF DAILY 771
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_771nm ~ mean_GOSAT_SIF_Daily_771nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Daily_771nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Daily_771nm, r_SIF_Daily_771nm, RMSE_SIF_Daily_771nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Daily 771nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Daily 771nm")), 1, 1.5, cex = 1)
box()

# SIF DAILY 740
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_740nm ~ mean_GOSAT_SIF_Daily_740nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Daily_740nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Daily_740nm, r_SIF_Daily_740nm, RMSE_SIF_Daily_740nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Daily 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Daily 740nm")), 1, 1.5, cex = 1)
box()

# SIF Relative 757
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_757nm ~ mean_GOSAT_SIF_Relative_757nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Relative_757nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Relative_757nm, r_SIF_Relative_757nm, RMSE_SIF_Relative_757nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Relative 757nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Relative 757nm")), 1, 1.5, cex = 1)
box()

# SIF Relative 771
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_771nm ~ mean_GOSAT_SIF_Relative_771nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Relative_771nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Relative_771nm, r_SIF_Relative_771nm, RMSE_SIF_Relative_771nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Relative 771nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Relative 771nm")), 1, 1.5, cex = 1)
box()

# dev.off()
#endregion

#region ############# Plotting GOSAT and OCO footprints ############

  # Transform OCO data to SpatialPolygonDataFrame
polydf_oco2 <- build_polyDF(veg_matched_oco2_sounding_list)
shapefile(polydf_oco2, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Veg_Matched_OCO2.shp", overwrite = TRUE)

coords_gosat   <- as.data.frame(cbind(veg_matched_gosat_sounding_list$longitude, veg_matched_gosat_sounding_list$latitude))
colnames(coords_gosat) <- c("longitude", "latitude")
point_df_gosat <- SpatialPointsDataFrame(coords_gosat, proj4string = CRS("+init=epsg:4326"), data = veg_matched_gosat_sounding_list, coords.nrs = c(8, 9))
polydf_gosat   <- buffer(point_df_gosat, width = 5000, dissolve = FALSE) # Create SpatialPolygonDataFrame with 5km radius buffer
shapefile(point_df_gosat, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Veg_Matched_GOSAT_points.shp", overwrite = TRUE)
shapefile(polydf_gosat, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Veg_Matched_GOSAT.shp", overwrite = TRUE)

coords_gosat   <- as.data.frame(cbind(n_vza_time_veg_matched_sounding_means$longitude, n_vza_time_veg_matched_sounding_means$latitude))
colnames(coords_gosat) <- c("longitude", "latitude")
point_df_gosat <- SpatialPointsDataFrame(coords_gosat, proj4string = CRS("+init=epsg:4326"), data = n_vza_time_veg_matched_sounding_means, coords.nrs = c(8, 9))
shapefile(point_df_gosat, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/All_Filtered_Matched_GOSAT_points.shp", overwrite = TRUE)


plot(polydf_oco2)
plot(polydf_gosat[1, ], add = TRUE)

#endregion


##### HISTOGRAM OF TIME ######

hist(veg_matched_sounding_means$Delta_Time, breaks = "quarters", freq = TRUE, format = "%Y-%m", right = FALSE)
