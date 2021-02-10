library(raster)
library(rgdal)
library(ncdf4)

options(scipen = 999)

files_gosat_2015 <- list.files(path = "C:/Russell/Projects/Geometry/Data/gosat/2015", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_gosat_2016 <- list.files(path = "C:/Russell/Projects/Geometry/Data/gosat/2016", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_gosat_2017 <- list.files(path = "C:/Russell/Projects/Geometry/Data/gosat/2017", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_gosat_2018 <- list.files(path = "C:/Russell/Projects/Geometry/Data/gosat/2018", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_gosat_2019 <- list.files(path = "C:/Russell/Projects/Geometry/Data/gosat/2019", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)

files_oco2_2015  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2/2015", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_oco2_2016  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2/2016", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_oco2_2017  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2/2017", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_oco2_2018  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2/2018", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_oco2_2019  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2/2019", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)

cosd <- function(degrees) {
  radians <- cos(degrees * pi / 180)
  return(radians)
}
sind <- function(degrees) {
  radians <- sin(degrees * pi / 180)
  return(radians)
}
compute_phase_angle <- function(sat) {
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
    temp        <- ncvar_get(nc, "Meteo/temperature_skin")
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
    df <- data.frame("SoundingID" = id, "MeasurementMode" = mode, "cloud_flag_abp" = cloud_flag, "Quality_Flag" = q_flag,
                    "Delta_Time" = as.POSIXct(time, origin = "1990-01-01", tz = "UTC"), "longitude" = lon_center, "latitude" = lat_center,
                    "SIF_Daily_740nm" = sif740_D, "SIF_Daily_757nm" = sif757_D, "SIF_Daily_771nm" = sif771_D,
                    "SIF_740nm" = sif740,
                    "SIF_757nm" = sif757, "SIF_Relative_757nm" = sif757_R, "continuum_radiance_757nm" = rad757,
                    "SIF_771nm" = sif771, "SIF_Relative_771nm" = sif771_R, "continuum_radiance_771nm" = rad771,
                    "IGBP_index" = igbp,
                    "SZA" = sza, "SAz" = saa, "VZA" = vza, "VAz" = vaa, "PA" = pa, "RAz" = raa)
  } else {
    df <- data.frame("SoundingID" = id, "MeasurementMode" = mode, "cloud_flag_abp" = cloud_flag,
                "Quality_Flag_P" = q_flag[1, ], "Quality_Flag_S" = q_flag[2, ],
                "Delta_Time" = as.POSIXct(time, origin = "1990-01-01", tz = "UTC"), "longitude" = lon_center, "latitude" = lat_center,
                "SIF_Daily_740nm_P" = sif740_D[1, ], "SIF_Daily_757nm_P" = sif757_D[1, ], "SIF_Daily_771nm_P" = sif771_D[1, ],
                "SIF_Daily_740nm_S" = sif740_D[2, ], "SIF_Daily_757nm_S" = sif757_D[2, ], "SIF_Daily_771nm_S" = sif771_D[2, ],
                "SIF_740nm_P" = sif740[1, ],
                "SIF_740nm_S" = sif740[2, ],
                "SIF_757nm_P" = sif757[1, ], "SIF_Relative_757nm_P" = sif757_R[1, ], "continuum_radiance_757nm_P" = rad757[1, ],
                "SIF_757nm_S" = sif757[2, ], "SIF_Relative_757nm_S" = sif757_R[2, ], "continuum_radiance_757nm_S" = rad757[2, ],
                "SIF_771nm_P" = sif771[1, ], "SIF_Relative_771nm_P" = sif771_R[1, ], "continuum_radiance_771nm_P" = rad771[1, ],
                "SIF_771nm_S" = sif771[2, ], "SIF_Relative_771nm_S" = sif771_R[2, ], "continuum_radiance_771nm_S" = rad771[2, ],
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


gosat_annual_2019 <- data.frame()

for (i in 1:length(files_gosat_2019)) {
    date <- substr(basename(files_gosat_2019[i]), 13, 18)
    print(paste0("Working on file number ", i, " and date ", date, "."))

    df <- build_data(files_gosat_2019[i], "GOSAT")
    df <- subset_flags(df, 0, 0, 0, "GOSAT")
    gosat_annual_2019 <- rbind(gosat_annual_2019, df)
}

gosat_annual_2019$SIF_Daily_740nm_Mean <- (gosat_annual_2019$SIF_Daily_740nm_P + gosat_annual_2019$SIF_Daily_740nm_S) / 2

coords                     <- as.data.frame(cbind(gosat_annual_2019$longitude, gosat_annual_2019$latitude))
colnames(coords)           <- c("longitude", "latitude")
point_df_gosat_annual_2019 <- SpatialPointsDataFrame(coords, proj4string = CRS("+init=epsg:4326"), data = subset(gosat_annual_2019, select = -c(longitude, latitude)), coords.nrs = c(7, 8))

r <- raster(ncols = 90, nrows = 45, xmn = -180, xmx = 180, ymn = -90, ymx = 90, crs = "+init=epsg:4326")
raster_gosat_annual_2019 <- rasterize(point_df_gosat_annual_2019, r, field = "SIF_Daily_740nm_Mean", fun = mean)
writeRaster(raster_gosat_annual_2019, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2019.tif", NAflag = -9999)


oco2_annual_2015 <- data.frame()

for (i in 1:length(files_oco2_2015)) {
    date <- substr(basename(files_oco2_2015[i]), 12, 17)
    print(paste0("Working on file number ", i, " and date ", date, "."))

    df <- build_data(files_oco2_2015[i], "OCO")
    df <- subset_flags(df, 0, 0, 0, "OCO")
    oco2_annual_2015 <- rbind(oco2_annual_2015, df)
}

coords                     <- as.data.frame(cbind(oco2_annual_2015$longitude, oco2_annual_2015$latitude))
colnames(coords)           <- c("longitude", "latitude")
point_df_oco2_annual_2015 <- SpatialPointsDataFrame(coords, proj4string = CRS("+init=epsg:4326"), data = subset(oco2_annual_2015, select = -c(longitude, latitude)), coords.nrs = c(6, 7))

r <- raster(ncols = 90, nrows = 45, xmn = -180, xmx = 180, ymn = -90, ymx = 90, crs = "+init=epsg:4326")
raster_oco2_annual_2015 <- rasterize(point_df_oco2_annual_2015, r, field = "SIF_Daily_740nm", fun = mean)
writeRaster(raster_oco2_annual_2015, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/OCO2_Annual_Mean_740Daily_2015.tif", NAflag = -9999)



