library(raster)
library(rgdal)
library(ncdf4)

options(scipen = 999)

files_gosat_june_2019 <- list.files(path = "C:/Russell/Projects/Geometry/Data/gosat/2019/06", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)

# files_oco2_2019  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2/2019", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_oco2_june_2020  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco2/2020/06", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)

# files_oco3_2019  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco3/2019", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)
files_oco3_june_2020  <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco3/2020/06", pattern = "*.nc4", full.names = TRUE, recursive = TRUE)

build_data <- function (input_file, platform) {
  env <- new.env()
  nc <- nc_open(input_file) # Open file
  # Metadata
  mode  <- ncvar_get(nc, "Metadata/MeasurementMode") # 0=Nadir, 1=Glint, 2=Target, 3=AreaMap, 4=Transition
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

  # Build Dataframe ####
  if (platform == "OCO"){
    df <- data.frame("MeasurementMode" = mode, "cloud_flag_abp" = cloud_flag, "Quality_Flag" = q_flag,
                    "longitude" = lon_center, "latitude" = lat_center,
                    "SIF_Daily_740nm" = sif740_D, "SIF_Daily_757nm" = sif757_D, "SIF_Daily_771nm" = sif771_D,
                    "SIF_740nm" = sif740,
                    "SIF_757nm" = sif757, "SIF_Relative_757nm" = sif757_R, "continuum_radiance_757nm" = rad757,
                    "SIF_771nm" = sif771, "SIF_Relative_771nm" = sif771_R, "continuum_radiance_771nm" = rad771)
  } else {
    df <- data.frame("MeasurementMode" = mode, "cloud_flag_abp" = cloud_flag,
                "Quality_Flag_P" = q_flag[1, ], "Quality_Flag_S" = q_flag[2, ],
                "longitude" = lon_center, "latitude" = lat_center,
                "SIF_Daily_740nm_P" = sif740_D[1, ], "SIF_Daily_757nm_P" = sif757_D[1, ], "SIF_Daily_771nm_P" = sif771_D[1, ],
                "SIF_Daily_740nm_S" = sif740_D[2, ], "SIF_Daily_757nm_S" = sif757_D[2, ], "SIF_Daily_771nm_S" = sif771_D[2, ],
                "SIF_740nm_P" = sif740[1, ],
                "SIF_740nm_S" = sif740[2, ],
                "SIF_757nm_P" = sif757[1, ], "SIF_Relative_757nm_P" = sif757_R[1, ], "continuum_radiance_757nm_P" = rad757[1, ],
                "SIF_757nm_S" = sif757[2, ], "SIF_Relative_757nm_S" = sif757_R[2, ], "continuum_radiance_757nm_S" = rad757[2, ],
                "SIF_771nm_P" = sif771[1, ], "SIF_Relative_771nm_P" = sif771_R[1, ], "continuum_radiance_771nm_P" = rad771[1, ],
                "SIF_771nm_S" = sif771[2, ], "SIF_Relative_771nm_S" = sif771_R[2, ], "continuum_radiance_771nm_S" = rad771[2, ])
  }
  df <- na.omit(df) # drop rows that contain an NA anywhere
  rownames(df) <- 1:nrow(df)

  nc_close(nc)
  return(df)
}
subset_flags <- function(df, mode, flag_cloud, flag_qc, platform) {
  if (platform == "OCO"){
    if (is.na(mode)) {
      # All modes, Do nothing
    } else if (mode == 0) {
      lab_mode <- "NADIR"
      df <- subset(df, MeasurementMode == 0)
    } else if (mode == 1) {
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
    if (is.na(flag_qc)) {
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
  if (platform == "GOSAT") {
      # Mode GOSAT
    if (is.na(mode)) {
      # All modes, Do nothing
    } else if (mode == 0) {
      lab_mode <- "OB1D"
      df <- subset(df, MeasurementMode == 0)
    } else if (mode == 1) {
      lab_mode <- "OB2D"
      df <- subset(df, MeasurementMode == 1)
    } else if (mode == 2) {
      lab_mode <- "SPOD"
      df <- subset(df, MeasurementMode == 2)
    }
        # QC GOSAT
    if (is.na(flag_qc)) {
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

  # list_return <- list("lab_mode" = lab_mode, "lab_cloud" = lab_cloud, "lab_qc" = lab_qc)
  # list2env(list_return, .GlobalEnv)
  return(df)
}


####### GOSAT ############

gosat_june_2019 <- data.frame()

for (i in 1:length(files_gosat_june_2019)) {
    date <- substr(basename(files_gosat_june_2019[i]), 13, 18)
    print(paste0("Working on file number ", i, " and date ", date, "."))

    df <- build_data(files_gosat_june_2019[i], "GOSAT")
    df <- subset_flags(df, NA, 0, 10, "GOSAT")
    gosat_june_2019 <- rbind(gosat_june_2019, df)
}

gosat_june_2019$SIF_Daily_757nm_Mean <- (gosat_june_2019$SIF_Daily_757nm_P + gosat_june_2019$SIF_Daily_757nm_S) / 2

coords                     <- as.data.frame(cbind(gosat_june_2019$longitude, gosat_june_2019$latitude))
colnames(coords)           <- c("longitude", "latitude")
point_df_gosat_june_2019 <- SpatialPointsDataFrame(coords, proj4string = CRS("+init=epsg:4326"), data = subset(gosat_june_2019, select = -c(longitude, latitude)), coords.nrs = c(7, 8))

r                        <- raster(ncols = 720, nrows = 360, xmn = -180, xmx = 180, ymn = -90, ymx = 90, crs = "+init=epsg:4326")
raster_gosat_june_2019 <- rasterize(point_df_gosat_june_2019, r, field = "SIF_Daily_757nm_Mean", fun = mean)

plot(raster_gosat_june_2019, zlim = c(0,0.8))

writeRaster(raster_gosat_june_2019, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/GOSAT_June_Mean_757Daily_2019_0.50.tif", NAflag = -9999, overwrite = TRUE)



####### OCO2 ############

r          <- raster(ncols = 720, nrows = 360, xmn = -180, xmx = 180, ymn = -90, ymx = 90, crs = "+init=epsg:4326")
oco2_stack <- stack()

for (i in 1:length(files_oco2_june_2020)) {
    date <- substr(basename(files_oco2_june_2020[i]), 12, 17)
    print(paste0("Working on file number ", i, " and date ", date, "."))

    df <- build_data(files_oco2_june_2020[i], "OCO")
    df <- subset_flags(df, NA, 0, 10, "OCO")
    
    if (nrow(df) != 0) {
      coords           <- as.data.frame(cbind(df$longitude, df$latitude))
      colnames(coords) <- c("longitude", "latitude")
      point_df         <- SpatialPointsDataFrame(coords, proj4string = CRS("+init=epsg:4326"), data = subset(df, select = -c(longitude, latitude)), coords.nrs = c(6, 7))
      
      oco2_raster <- rasterize(point_df, r, field = "SIF_Daily_757nm", fun = mean)
      oco2_stack  <- stack(oco2_raster, oco2_stack)
    }
}

oco2_stack_mean <- mean(oco2_stack, na.rm = TRUE)

plot(oco2_stack_mean, zlim = c(0,0.8))

writeRaster(oco2_stack_mean, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/OCO2_June_Mean_757Daily_2020_0.50.tif", NAflag = -9999, overwrite = TRUE)


####### OCO3 ############

r          <- raster(ncols = 720, nrows = 360, xmn = -180, xmx = 180, ymn = -90, ymx = 90, crs = "+init=epsg:4326")
oco3_stack <- stack()

for (i in 1:length(files_oco3_june_2020)) {
  date <- substr(basename(files_oco3_june_2020[i]), 12, 17)
  print(paste0("Working on file number ", i, " and date ", date, "."))
  
  df <- build_data(files_oco3_june_2020[i], "OCO")
  df <- subset_flags(df, NA, 0, 10, "OCO")
  
  if (nrow(df) != 0) {
    coords           <- as.data.frame(cbind(df$longitude, df$latitude))
    colnames(coords) <- c("longitude", "latitude")
    point_df         <- SpatialPointsDataFrame(coords, proj4string = CRS("+init=epsg:4326"), data = subset(df, select = -c(longitude, latitude)), coords.nrs = c(6, 7))
    
    oco3_raster <- rasterize(point_df, r, field = "SIF_Daily_757nm", fun = mean)
    oco3_stack  <- stack(oco3_raster, oco3_stack)
  }
}

oco3_stack_mean <- mean(oco3_stack, na.rm = TRUE)

plot(oco3_stack_mean, zlim = c(0,0.8))

writeRaster(oco3_stack_mean, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/OCO3_June_Mean_757Daily_2020_0.50.tif", NAflag = -9999, overwrite = TRUE)



