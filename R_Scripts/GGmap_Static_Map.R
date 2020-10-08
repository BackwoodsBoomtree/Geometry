library(ncdf4)
library(viridis)
library(ggmap)
library(ggplot2)
library(grid)
library(sp)
library(gtools)
library(raster)
library(rgdal)


options(scipen = 999)
options(digits = 14)

input_dir <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco3/niwot", full.names = TRUE, pattern = "*.nc4")
# input_dir <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco3/ecostress_us_syv", full.names = TRUE, pattern = "*.nc4")
# input_dir <- list.files(path = "C:/Russell/Projects/Geometry/Data/oco3/ATTO_incorrect", full.names = TRUE, pattern = "*.nc4")
# input_dir <- list.files(path = "C:/Russell/Projects/Geometry/oco3/Lamont", full.names = TRUE, pattern = "*.nc4")
target_list <- read.csv("C:/Russell/Projects/Geometry/Data/oco3/site_list/oco3_targets.csv")

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
find_site <- function(target_list) {
  # Find target
  target_list[, 3] <- gsub("POINT\\(", "", target_list[, 3])
  target_list[, 3] <- gsub("\\)", "", target_list[, 3])
  target_list_locs <- t(data.frame(strsplit(target_list[, 3], " ")))
  target_list$x <- as.numeric(target_list_locs[, 1])
  target_list$y <- as.numeric(target_list_locs[, 2])
  target_list[, 3:4] <- NULL
  colnames(target_list) <- c("Target_ID", "Target_Name", "lon", "lat")
  remove(target_list_locs)
  coordinates(target_list) <- ~ lon + lat
  # Set the projection of the SpatialPointsDataFrame using the projection of the shapefile
  proj4string(target_list) <- proj4string(poly_df)
  test <- over(target_list, poly_df)
  test2 <- over(poly_df, target_list)
}
build_data <- function (input_file) {
  env <- new.env()
  df <- nc_open(input_file) # Open file
  # Metadata
  mode <- ncvar_get(df, "Metadata/MeasurementMode") # 0=Nadir, 1=Glint, 2=Target, 3=AreaMap, 4=Transition
  orbit <- ncvar_get(df, "Metadata/OrbitId")
  id <- ncvar_get(df, "Metadata/SoundingId") # "YYYYMMDDHHMMSS"
  time <- ncvar_get(df, "Delta_Time") # seconds since 1990-01-01 00:00:00 UTC
  igbp <- ncvar_get(df, "Science/IGBP_index")
  percent_cover <- ncvar_get(df, "Science/sounding_land_fraction")
  # SIF and Flag Data
  cloud_flag <- ncvar_get(df, "Cloud/cloud_flag_abp") # 0 - \"Classified clear\", 1 - \"Classified cloudy\", 2 - \"Not classified\", all other values undefined; not used in SIF processing
  q_flag <- ncvar_get(df, "Quality_Flag") # 0 = best (passes quality control + cloud fraction = 0.0); 1 = good (passes quality control); 2 = bad (failed quality control); -1 = not investigated
  sif740_D <- ncvar_get(df, "Daily_SIF_740nm")
  sif757_D <- ncvar_get(df, "Daily_SIF_757nm")
  sif771_D <- ncvar_get(df, "Daily_SIF_771nm")
  sif740 <- ncvar_get(df, "SIF_740nm")
  sif757 <- ncvar_get(df, "Science/SIF_757nm")
  sif771 <- ncvar_get(df, "Science/SIF_771nm")
  sif740_U <- ncvar_get(df, "SIF_Uncertainty_740nm")
  sif757_U <- ncvar_get(df, "Science/SIF_Uncertainty_757nm")
  sif771_U <- ncvar_get(df, "Science/SIF_Uncertainty_771nm")
  # Radiance
  rad757 <- ncvar_get(df, "Science/continuum_radiance_757nm")
  rad771 <- ncvar_get(df, "Science/continuum_radiance_771nm")
  # Geo Data
  lon_corners <- ncvar_get(df, "Geolocation/footprint_longitude_vertices")
  lat_corners <- ncvar_get(df, "Geolocation/footprint_latitude_vertices")
  sza <- ncvar_get(df, "SZA")
  saa <- ncvar_get(df, "SAz")
  vza <- ncvar_get(df, "VZA")
  vaa <- ncvar_get(df, "VAz")
  # Relative azimuth angle
  raa <- abs(saa - vaa)
  for (i in 1:length(raa)) {
    if (raa[i] > 180) {
      raa[i] <- abs(raa[i] - 360)
    }
  }
  # Meteo
  humidity <- ncvar_get(df, "Meteo/specific_humidity")
  surface_pressure <- ncvar_get(df, "Meteo/surface_pressure")
  temp_skin <- ncvar_get(df, "Meteo/temperature_skin")
  temp_2m <- ncvar_get(df, "Meteo/temperature_two_meter")
  vpd <- ncvar_get(df, "Meteo/vapor_pressure_deficit")
  wind <- ncvar_get(df, "Meteo/wind_speed")
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

  df <- data.frame("sid" = id, "mode" = mode, "orbit" = orbit, "cloud_flag" = cloud_flag, "quality_flag" = q_flag,
                   "time" = as.POSIXct(time, origin = "1990-01-01", tz = "UTC"),
                   "sif740_D" = sif740_D, "sif757_D" = sif757_D, "sif771_D" = sif771_D,
                   "sif740" = sif740, "sif740_U" = sif740_U,
                   "sif757" = sif757, "sif757_U" = sif757_U, "rad757" = rad757,
                   "sif771" = sif771, "sif771_U" = sif771_U, "rad771" = rad771,
                   "lon_1" = lon1, "lon_2" = lon2, "lon_3" = lon3, "lon_4" = lon4,
                   "lat_1" = lat1, "lat_2" = lat2, "lat_3" = lat3, "lat_4" = lat4,
                   "sza" = sza, "saa" = saa, "vza" = vza, "vaa" = vaa, "pa" = pa, "raa" = raa,
                   "humidity" = humidity, "surface_pressure" = surface_pressure, "temp_skin" = temp_skin,
                   "temp_2m" = temp_2m, "vpd" = vpd, "wind" = wind, "igbp_class" = igbp, "percent_cover" = percent_cover)
  df <- na.omit(df) # drop rows that contain an NA anywhere
  # Get time for plot labels and filenaming
  lab_time <<- as.POSIXct((as.numeric(max(df$time)) + as.numeric(min(df$time))) / 2, origin = "1970-01-01", tz = "UTC")
  return(df)
}
build_data_multi <- function (input_dir) {
  for (f in 1:(length(input_dir))) {
    df_temp <- build_data(input_dir[f])
    df_temp <- subset_flags(df_temp, mode, flag_cloud, flag_qc)
    df_temp <- subset_location(df_temp, lat_min, lat_max, lon_min, lon_max)
    df_temp <- subset_cover(df_temp, igbp)
    if (f == 1) {
      df <- df_temp
      assign(paste0("df", f), df, envir = .GlobalEnv) # Return each df individually
    } else {
      df <- rbind(df, df_temp)
      assign(paste0("df", f), df_temp, envir = .GlobalEnv) # Return each df individually
    }
  }
  return(df)
}
subset_flags <- function(df, mode, flag_cloud, flag_qc) {
  # mode
  if (mode == 0) {
    lab_mode <- "NADIR"
    df <- subset(df, mode == 0)
  }else if (mode == 1) {
    lab_mode <- "GLINT"
    df <- subset(df, mode == 1)
  } else if (mode == 2) {
    lab_mode <- "TARGET"
    df <- subset(df, mode == 2)
  } else if (mode == 3) {
    lab_mode <- "SAM"
    df <- subset(df, mode == 3)
  } else if (mode == 4) {
    lab_mode <- "TRANSITION"
    df <- subset(df, mode == 4)
  } else if (mode == 5) {
    lab_mode <- "SAM & Target"
    df <- subset(df, mode == 2 | mode == 3)
  }
  # Cloud
  if (is.na(flag_cloud)) {
    lab_cloud <- "None"
  } else if (flag_cloud == 0) {
    lab_cloud <- "Clear"
    df <- subset(df, cloud_flag == 0)
  } else if (flag_cloud == 1) {
    lab_cloud <- "Cloudy"
    df <- subset(df, cloud_flag == 1)
  } else if (flag_cloud == 2) {
    lab_cloud <- "Not_Classified"
    df <- subset(df, cloud_flag == 2)
  } else if (flag_cloud == 10) {
    lab_cloud <- "Clear?Cloudy"
    df <- subset(df, cloud_flag <= 1)
  }
  # QC
  if (is.na(flag_qc)){
    lab_qc <- "None"
  } else if (flag_qc == 0) {
    lab_qc <- "Best"
    df <- subset(df, quality_flag == 0)
  } else if (flag_qc == 1) {
    lab_qc <- "Good"
    df <- subset(df, quality_flag == 1)
  } else if (flag_qc == 2) {
    lab_qc <- "Bad"
    df <- subset(df, quality_flag == 2)
  } else if (flag_qc == 10) {
    lab_qc <- "Best/Good"
    df <- subset(df, quality_flag <= 1)
  } else if (flag_qc == -1) {
    lab_qc <- "Not_Investigated"
    df <- subset(df, quality_flag == -1)
  }
  lab_class <- "All" # replaced by subset_cover function if it is used to subset land cover
  list_return <- list("lab_mode" = lab_mode, "lab_cloud" = lab_cloud, "lab_qc" = lab_qc, "lab_class" = lab_class)
  list2env(list_return, .GlobalEnv)
  return(df)
}
subset_location <- function (df, lat_min, lat_max, lon_min, lon_max) {
  df <- subset(df, lat_1 > lat_min & lat_1 < lat_max & lon_1 > lon_min & lon_1 < lon_max)
  return(df)
}
subset_orbit <- function (df, orbit_num) {
  df <- subset(df, orbit == orbit_num)
  return(df)
}
subset_cover <- function (df, igbp, percent) {
  if (!is.na(percent)) {
    df <- subset(df, percent_cover >= percent) # Filter by percent
    lab_percent <<- paste0(percent, "%")
  } else if (is.na(percent)) {
    lab_percent <<- paste0("")
  }
  if (is.na(igbp)) {
    lab_class <- "All"
  } else if (igbp == 1) {
    lab_class <- "Evergreen NF"
    df <- subset(df, igbp_class == 1)
  } else if (igbp == 2) {
    lab_class <- "Evergreen BF"
    df <- subset(df, igbp_class == 2)
  } else if (igbp == 3) {
    lab_class <- "Deciduous NF"
    df <- subset(df, igbp_class == 3)
  } else if (igbp == 4) {
    lab_class <- "Deciduous BF"
    df <- subset(df, igbp_class == 4)
  } else if (igbp == 5) {
    lab_class <- "Mixed Forest"
    df <- subset(df, igbp_class == 5)
  } else if (igbp == 6) {
    lab_class <- "Closed Shrub"
    df <- subset(df, igbp_class == 6)
  } else if (igbp == 7) {
    lab_class <- "Open Shrub"
    df <- subset(df, igbp_class == 7)
  } else if (igbp == 8) {
    lab_class <- "Woody Savanna"
    df <- subset(df, igbp_class == 8)
  } else if (igbp == 9) {
    lab_class <- "Savanna"
    df <- subset(df, igbp_class == 9)
  } else if (igbp == 10) {
    lab_class <- "Grassland"
    df <- subset(df, igbp_class == 10)
  } else if (igbp == 11) {
    lab_class <- "Perm Wetland"
    df <- subset(df, igbp_class == 11)
  } else if (igbp == 12) {
    lab_class <- "Cropland"
    df <- subset(df, igbp_class == 12)
  } else if (igbp == 13) {
    lab_class <- "Urban/Built-up"
    df <- subset(df, igbp_class == 13)
  } else if (igbp == 14) {
    lab_class <- "Crop/Veg Mosaic"
    df <- subset(df, igbp_class == 14)
  } else if (igbp == 15) {
    lab_class <- "Snow/Ice"
    df <- subset(df, igbp_class == 15)
  } else if (igbp == 16) {
    lab_class <- "Barren"
    df <- subset(df, igbp_class == 16)
  } else if (igbp == 17) {
    lab_class <- "Water"
  }
  list_return <- list("lab_class" = lab_class)
  list2env(list_return, .GlobalEnv)
  return(df)
}
build_polyDF <- function(df) {
  #### Create Polygons for each row ####
  for (i in 1:nrow(df)) {
    x <- c(df$lon_1[i], df$lon_2[i], df$lon_3[i], df$lon_4[i], df$lon_1[i])
    y <- c(df$lat_1[i], df$lat_2[i], df$lat_3[i], df$lat_4[i], df$lat_1[i])
    #coords_df <- data.frame(cbind(x,y))
    #hull <- chull(coords_df)
    #hull <- c(hull, hull[1])
    #coords_df <- coords_df[match(hull, row.names(coords_df)),]
    #polly <- Polygon(coords_df)
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
build_laiCop <- function (cop_data, poly_df){
  # Input LAI raster is from Copernicus (1 km)
  lai_cop <- brick(cop_data, varname = "LAI")
  #extent(lai_cop) <- extent(-180, 180, -60, 80)
  # Sample LAI raster using OCO footprint polygons and calculate area weighted mean
  lai_cop <- extract(lai_cop, poly_df, weights = TRUE, fun = mean, na.rm = TRUE)
  lai_cop <- round(lai_cop, digits = 3)
  poly_df$lai_cop <- as.vector(lai_cop) # add LAI to shapefile
  return(poly_df)
}
build_laiERA5 <- function (era5_data, poly_df){
  # To get LAI for a pixel, we need to calculate from high and low vegetation cover
  high_cover <- brick(era5_data, varname = "cvh")
  low_cover <- brick(era5_data, varname = "cvl")
  lai_high <- brick(era5_data, varname = "lai_hv")
  lai_low <- brick(era5_data, varname = "lai_lv")
  lai_era5 <- (high_cover * lai_high) + (low_cover * lai_low)
  # Sample LAI raster using OCO footprint polygons and calculate area weighted mean
  lai_era5 <- round(as.vector(extract(lai_era5, poly_df, weights = TRUE, fun = mean, na.rm = TRUE)), digits = 3)
  poly_df$lai_era5 <- lai_era5 # add LAI to shapefile
  return(poly_df)
}
build_incoming_sw_ERA5 <- function (era5_data, poly_df){
  incoming_sw_era5 <- brick(era5_data, varname = "msdwswrfcs") # Mean surface downward short-wave radiation flux, clear sky
  incoming_direct_era5 <- brick(era5_data, varname = "msdrswrfcs") # Mean surface downward short-wave radiation flux, clear sky
  # Sample LAI raster using OCO footprint polygons and calculate area weighted mean
  incoming_sw_era5 <- round(as.vector(extract(incoming_sw_era5, poly_df, weights = TRUE, fun = mean, na.rm = TRUE)), digits = 3)
  incoming_direct_era5 <- round(as.vector(extract(incoming_direct_era5, poly_df, weights = TRUE, fun = mean, na.rm = TRUE)), digits = 3)
  incoming_diffuse_era5 <- incoming_sw_era5 - incoming_direct_era5
  # Add to shapefile
  poly_df$incoming_sw_era5 <- incoming_sw_era5
  poly_df$incoming_direct_era5 <- incoming_direct_era5
  poly_df$incoming_diffuse_era5 <- incoming_diffuse_era5
  return(poly_df)
}
build_evi <- function (evi_raster, poly_df){
  # Input EVI raster is from VPM input
  evi_raster <- raster(evi_raster)
  evi_raster <- projectRaster(evi_raster, crs = "+init=epsg:4326")
  # Sample LAI raster using OCO footprint polygons and calculate area weighted mean
  evi_awm <- extract(evi_raster, poly_df, weights = TRUE, fun = mean, na.rm = TRUE)
  evi_awm <- round((evi_awm / 10000), digits = 3)
  poly_df$evi <- as.vector(evi_awm) # add LAI to shapefile
  return(poly_df)
}
plot_data <- function (df, variable, save, site_name, output_dir, offset) {
  register_google(key = "AIzaSyDPeI_hkrch7DqhKmhFJKeADWBpAKJL3h4")

  # Labels for plotting
  lab_loc <- gsub("_", " ", site_name)
  lab_var <- sapply(gsub("_", " ", variable), toupper)
  if (variable == "igbp_class") {
    df$categorical <- df[[variable]]
    df$categorical <- gsub("\\<1\\>", "1 - Evergreen NF", df$categorical)
    df$categorical <- gsub("\\<2\\>", "2 - Evergreen BF", df$categorical)
    df$categorical <- gsub("\\<3\\>", "3 - Deciduous NF", df$categorical)
    df$categorical <- gsub("\\<4\\>", "4 - Deciduous BF", df$categorical)
    df$categorical <- gsub("\\<5\\>", "5 - Mixed Forest", df$categorical)
    df$categorical <- gsub("\\<6\\>", "6 - Closed Shrub", df$categorical)
    df$categorical <- gsub("\\<7\\>", "7 - Open Shrub", df$categorical)
    df$categorical <- gsub("\\<8\\>", "8 - Woody Savanna", df$categorical)
    df$categorical <- gsub("\\<9\\>", "9 - Savanna", df$categorical)
    df$categorical <- gsub("\\<10\\>", "10 - Grassland", df$categorical)
    df$categorical <- gsub("\\<11\\>", "11 - Perm Wetlands", df$categorical)
    df$categorical <- gsub("\\<12\\>", "12 - Cropland", df$categorical)
    df$categorical <- gsub("\\<13\\>", "13 - Urban/Built-up", df$categorical)
    df$categorical <- gsub("\\<14\\>", "14 - Crop/Veg Mosaic", df$categorical)
    df$categorical <- gsub("\\<15\\>", "15 - Snow/Ice", df$categorical)
    df$categorical <- gsub("\\<16\\>", "16 - Barren", df$categorical)
    df$categorical <- gsub("\\<17\\>", "17 - Water", df$categorical)
    cat_labels <- mixedsort(unique(df$categorical))
  } else if (variable == "cloud_flag") {
    df$categorical <- df[[variable]]
    df$categorical <- gsub("\\<0\\>", "0 - Clear", df$categorical)
    df$categorical <- gsub("\\<1\\>", "1 - Cloudy", df$categorical)
    df$categorical <- gsub("\\<2\\>", "2 - Not Classified", df$categorical)
    cat_labels <- mixedsort(unique(df$categorical))
  } else if (variable == "quality_flag") {
    df$categorical <- df[[variable]]
    df$categorical <- gsub("\\<0\\>", "0 - Best", df$categorical)
    df$categorical <- gsub("\\<1\\>", "1 - Good", df$categorical)
    df$categorical <- gsub("\\<2\\>", "2 - Bad", df$categorical)
    df$categorical <- gsub("\\<-1\\>", "-1 - Not Investigated", df$categorical)
    cat_labels <- mixedsort(unique(df$categorical))
  } else {
    df$categorical <- cut(df[[variable]], 7, include.lowest = TRUE, dig.lab = 3)
    cat_labels <- gsub("\\[|\\]|\\(|\\)", "", sort(unique(df$categorical))) # remove brackets
    cat_labels <- gsub("\\,", " - ", cat_labels)
  }
  # ggplot2 doesn't process SpatialPolygonDataFrame directly, so we must use fortify
  df$id <- rownames(df@data)
  df_fort <- fortify(df, region = "id") # id column here is a unique identifier for each row
  df_plot <- merge(df_fort, df@data, by = "id")
  #### Plot ####
  # Get min and max lat and lon
  x_min <- min(c(df$lon_1, df$lon_2, df$lon_3, df$lon_4))
  x_max <- max(c(df$lon_1, df$lon_2, df$lon_3, df$lon_4))
  y_min <- min(c(df$lat_1, df$lat_2, df$lat_3, df$lat_4))
  y_max <- max(c(df$lat_1, df$lat_2, df$lat_3, df$lat_4))
  # Make bounding box for calculating autozoom
  map_bb <- make_bbox(c(x_min, x_max), c(y_min, y_max))
  autozoom <- calc_zoom(map_bb) + offset
  # Center Point
  x_length <- abs(x_max - x_min)
  y_length <- abs(y_max - y_min)
  center_point <- c(lon = (x_min + (0.5 * x_length)), lat = (y_min + (0.5 * y_length)))
  # Main map
  map_main <- ggmap(get_map(location = c(lon = as.numeric(center_point[1]), lat = as.numeric(center_point[2])), source = "google", maptype="satellite", zoom = autozoom)) +
    geom_polygon(data = df_plot, aes(x = long, y = lat, group = group, fill = factor(categorical))) +
    scale_fill_manual(
      values = rev(plasma(length(unique(df$categorical)))),
      breaks = rev(mixedsort(unique(df$categorical))),
      labels = rev(cat_labels)) +
    labs(title = paste0(lab_loc, "\n", lab_time, " | Orbit ", df$orbit[1], " | Mode: ", lab_mode,
                        "\nCover: ", lab_class, " ", lab_percent, " | QC Filter: ", gsub("_", " ", lab_qc), " | Cloud Filter: ", gsub("_", " ", lab_cloud)),
                        fill = lab_var) +
    theme(axis.ticks = element_blank(), axis.title = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
          plot.title = element_text(hjust = 0.5), legend.justification = c(0, 1), legend.position = c(1.025, 1), legend.key.size = unit(1.5, 'lines'))
  # Inset map
  map_loc <- ggmap(get_map(location = c(lon = as.numeric(center_point[1]), lat = as.numeric(center_point[2])), source = "google", maptype = "satellite", zoom = 3)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    geom_point(aes(x = as.numeric(center_point[1]), y = as.numeric(center_point[2])), color = "white", size = 2, shape = 21, fill = "#c51b7d", show.legend = FALSE) +
    theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), plot.title = element_blank(), panel.border = element_rect(colour = "black", fill = NA))
  # Histogram or Bar
  if (variable == "igbp_class"){
    h <- ggplot(df@data, aes(x = factor(igbp_class), fill=factor(igbp_class))) + xlab(lab_var) +
      geom_bar(show.legend = FALSE) +
      scale_fill_manual(
        values = plasma(length(unique(df$categorical)))) +
      scale_y_continuous(expand = expansion(mult = c(0, .1))) +
      theme(panel.border = element_rect(colour = "black", fill=NA), panel.background = element_blank(),
            axis.title.y = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0, 0, 0, -0.2), "cm"),
            axis.text.x = element_text(margin=unit(c(0, 0, 0.1, 0), "cm")),
            axis.text.y = element_blank(),
            axis.title.x = element_text(margin=unit(c(0, 0, 0, 0), "cm")))
  } else if (variable == "cloud_flag"){
    h <- ggplot(df@data, aes(x = factor(cloud_flag), fill=factor(cloud_flag))) + xlab(lab_var) +
      geom_bar(show.legend = FALSE) +
      scale_fill_manual(
        values = plasma(length(unique(df$categorical)))) +
      scale_y_continuous(expand = expansion(mult = c(0, .1))) +
      theme(panel.border = element_rect(colour = "black", fill = NA), panel.background = element_blank(),
            axis.title.y = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0, 0, 0, -0.2), "cm"),
            axis.text.x = element_text(margin = unit(c(0, 0, 0.1, 0), "cm")),
            axis.text.y = element_blank(),
            axis.title.x = element_text(margin = unit(c(0, 0, 0, 0), "cm")))
  } else if (variable == "quality_flag"){
    h <- ggplot(df@data, aes(x=factor(quality_flag), fill=factor(quality_flag))) + xlab(lab_var) +
      geom_bar(show.legend = FALSE) +
      scale_fill_manual(
        values = plasma(length(unique(df$categorical)))) +
      scale_y_continuous(expand = expansion(mult = c(0, .1))) +
      theme(panel.border = element_rect(colour = "black", fill = NA), panel.background = element_blank(),
            axis.title.y = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0, 0, 0, -0.2), "cm"),
            axis.text.x = element_text(margin = unit(c(0, 0, 0.1, 0), "cm")),
            axis.text.y = element_blank(),
            axis.title.x = element_text(margin = unit(c(0, 0, 0, 0), "cm")))
  } else {
    h <- ggplot(df@data, aes_string(x = variable)) + xlab(lab_var) +
      geom_histogram(aes(y = ..density..), binwidth = 0.25, color = "black", fill = "gray85") +
      geom_density(alpha = .2, fill = "#FF6666") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0, .1))) +
      theme(panel.border = element_rect(colour = "black", fill = NA), panel.background = element_blank(),
            axis.title.y = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0, 0, 0, -0.2), "cm"),
            axis.text.x = element_text(margin = unit(c(0, 0, 0.1, 0), "cm")),
            axis.text.y = element_blank(),
            axis.title.x = element_text(margin = unit(c(0, 0, 0, 0), "cm")))
  }
  # Viewports
  vp_main <- viewport(x = 0.5, y = 0.5)
  vp_map <- viewport(x = 0.375, y = 0.5, width = 1)
  vp1 <- viewport(width = 0.19, height = 0.15, x = 0.87, y = 0.375)
  vp2 <- viewport(width = 0.23, height = 0.23, x = 0.87, y = 0.15)
  if (save) {
    pdf(paste0(output_dir, site_name, "_", variable, "_", lab_mode, "_", df$orbit[1], "_Cloud", gsub("/", "", lab_cloud), "_QC",
               gsub("/", "", lab_qc), "_Cover", gsub(" ", "", lab_class), ".pdf"), width = 7.5, height = 6.25, compress = FALSE)
    print(map_main, vp = vp_map)
    print(h, vp = vp1)
    print(map_loc, vp = vp2)
    dev.off()
  }
  print(map_main, vp = vp_map)
  print(h, vp = vp1)
  print(map_loc, vp = vp2)
}

#### PLOTTING SITES ####
df <- build_data(input_dir[1])

# Args: input df, mode, cloud flag, qc flag
# mode: 0 = Nadir; 1 = Glint; 2 = Target; 3 = SAM; 4 = Transition; 5 = SAM & Target
# cloud flag: NA = no filter; 0 = clear; 1 = cloudy; 2 = Not classified; 10 = clear and cloudy
# qc flag: NA = no filter; 0 = best; 1 = good; 2 = bad; -1 = not investigated; 10 = best and good
df <- subset_flags(df, 3, 0, 0)

# min lat, max lat, min lon, max lon
df <- subset_location(df, 39, 42, -108, -104) # niwot
# df <- subset_location(df, 0, 5, -61, -57) # ATTO - incorrect
# f <- subset_location(df, 35, 37, 139, 141) # Tokyo
# df <- subset_location(df, 34, 38, -100, -94) # Lamont

# IGBP number and percent
df <- subset_cover(df, 1, NA) # niwot
# df <- subset_cover(df, NA, NA) # ecostress_us_syv
# df <- subset_cover(df, 2, 100) # Amazon
# df <- subset_cover(df, "None") # Lamont

# Orbit number
# df_6283 <- subset_orbit(df, 6283)
df_6287 <- subset_orbit(df, 6287)

poly_df <- build_polyDF(df_6283) # Build shapefile

# Add LAI to shapefile
poly_df <- build_laiCop("C:/Russell/Projects/Geometry/Data/lai/c_gls_LAI-RT0_202006300000_GLOBE_PROBAV_V2.0.1.nc", poly_df)

# Add Shortwave Radiation Downward at the surface
poly_df <- build_incoming_sw_ERA5("C:/Russell/Projects/Geometry/Data/era5/ERA5_SWDOWN_2020-06-12_1700.nc", poly_df) # niwot
# poly_df <- build_incoming_sw_ERA5("C:/Russell/Projects/Geometry/Data/era5/ERA5_SWDOWN_2020-06-17.nc", poly_df) # ecostress_us_syv
# poly_df <- build_incoming_sw_ERA5("C:/Russell/Projects/Geometry/Data/era5/ERA5_SWDOWN_2020-06-26.nc", poly_df) # ATTO - incorrect

# Add EVI to shapefile
# poly_df <- build_evi("C:/Russell/Projects/Geometry/Data/evi/GPP.2019177.h12v08.tif", poly_df)

# Arg: SpatialPolygonDF, variable of interest, save to file?, site name, output directory name, offset autozoom
plot_data(poly_df, "sif740", TRUE, "Niwot", "C:/Russell/Projects/Geometry/R_Scripts/Figures/", -1)
# plot_data(poly_df, "sif740", TRUE, "ecostress_us_syv", "C:/Russell/Projects/Geometry/R_Scripts/Figures/")
# plot_data(poly_df, "sif740", TRUE, "sif_ATTO_Tower_Manaus_Brazil_(incorrect)", "C:/Russell/Projects/Geometry/R_Scripts/Figures/")
# plot_data(poly_df, "sif740_D", FALSE, "val_tsukubaJp", "C:/Russell/R_Scripts/Geometry/")
# plot_data(poly_df, "sif740", TRUE, "val_lamontOK", "C:/Russell/R_Scripts/Geometry/")

# Export df to csv
write.csv(as.data.frame(poly_df), paste0("C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_niwot_", format(lab_time, format = '%Y-%m-%d', usetz = FALSE), "_", poly_df$orbit[1], ".csv"), row.names = FALSE)
# write.csv(as.data.frame(poly_df), paste0("C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_ecostress_us_syv_", format(lab_time, format = '%Y-%m-%d', usetz = FALSE), ".csv"), row.names = FALSE)
# write.csv(as.data.frame(poly_df), paste0("C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_ATTO_Tower_Manaus_Brazil_(incorrect)_", format(lab_time, format = '%Y-%m-%d', usetz = FALSE), ".csv"), row.names = FALSE)

# Export to shapefile
shapefile(poly_df, paste0("C:/Russell/Projects/Geometry/R_Scripts/SHP/sif_ATTO_Tower_Manaus_Brazil_(incorrect)_", format(lab_time, format = '%Y-%m-%d', usetz = FALSE), "_", poly_df$orbit[1], ".shp"), overwrite = TRUE)





#### PLOTTING MULTISITE ####

input_dir <- list.files(path = "C:/Russell/Projects/Geometry/oco3/ATTO_incorrect", full.names = TRUE, pattern = "*.nc4")
mode <- 3
flag_cloud <- 2
flag_qc <- NA
lat_min <- 0
lat_max <- 5
lon_min <- -62
lon_max <- -56
igbp <- 2
variable <- "sif740"
save <- TRUE
site_name <- "sif_ATTO_Tower_Manaus_Brazil_(incorrect)"
output_dir <- "C:/Russell/Projects/Geometry/"

df <- build_data_multi(input_dir)
poly_multi <- build_polyDF(df)
plot_data(poly_multi, variable, save, site_name, output_dir)

#### Multi vza ####
cols <- plasma(length(input_dir) * 3)
cols <- cols[c(TRUE, FALSE, FALSE)] # Remove every other color
par(xpd = TRUE)
plot(NA, xlim = c(0, 60), ylim = c(-2.75, 6), xlab = "vza", ylab = "sif740")
for (i in 1:length(input_dir)) {
  df_temp <- get(paste0("df", i))
  if (nrow(df_temp > 0)) { # Only execute if the df has data
    points(df_temp$vza, df_temp$sif740, col = cols[i])
    # polynomial line
    model <- lm(df_temp$sif740 ~ df_temp$vza + I(df_temp$vza ^ 2))
    myPredict <- predict(model)
    ix <- sort(df_temp$vza, index.return = TRUE)$ix
    lines(df_temp$vza[ix], myPredict[ix], col = cols[i], lwd = 4)
    ## rounded coefficients for better output
    cf <- round(coef(model), 4)
    eq <- paste0("y = ", round(cf[1], 2),
                ifelse(sign(cf[2]) == 1, " + ", " - "), abs(round(cf[2], 2)), "x ",
                ifelse(sign(cf[3]) == 1, " + ", " - "), format(abs(cf[3]), scientific=FALSE), "x^2")
    if (exists("legend_text") == FALSE) {
      legend_text <- paste0(eq, " | sza: ", round(min(df_temp$sza), digits = 1), " - ", round(max(df_temp$sza), digits = 1),
                            " | ", format(as.POSIXct(df_temp$time[1], format = "%Y/%m/%d %H:%M:%S"), format = "%Y/%m/%d"))
      text_cols <- cols[i]
    } else {
      legend_text <- c(legend_text, paste0(eq, " | sza: ", round(min(df_temp$sza), digits = 1), " - ", round(max(df_temp$sza), digits = 1),
                                          " | ", format(as.POSIXct(df_temp$time[1], format = "%Y/%m/%d %H:%M:%S"), format = "%Y/%m/%d")))
      text_cols <- c(text_cols, cols[i])
    }
  }
}
legend(0, 7.25, legend = legend_text, text.col = text_cols, text.font = 2)
remove(legend_text)

#### Multi pa ####
cols <- plasma(length(input_dir) * 3)
cols <- cols[c(TRUE, FALSE, FALSE)] # Remove every other color
plot(NA, xlim = c(0, 70), ylim = c(-2.75, 6), xlab = "pa", ylab = "sif740")
for (i in 1:length(input_dir)) {
  df_temp <- get(paste0("df", i))
  if (nrow(df_temp > 0)) { # Only execute if the df has data
    points(abs(df_temp$pa), df_temp$sif740, col = cols[i])
    # polynomial line
    model <- lm(df_temp$sif740 ~ abs(df_temp$pa) + I(abs(df_temp$pa) ^ 2))
    myPredict <- predict(model)
    ix <- sort(abs(df_temp$pa), index.return = TRUE)$ix
    lines(abs(df_temp$pa)[ix], myPredict[ix], col = cols[i], lwd = 4)
    ## rounded coefficients for better output
    cf <- round(coef(model), 4)
    eq <- paste0("y = ", round(cf[1], 2),
                ifelse(sign(cf[2]) == 1, " + ", " - "), abs(round(cf[2], 2)), "x ",
                ifelse(sign(cf[3]) == 1, " + ", " - "), format(abs(cf[3]), scientific = FALSE), "x^2")
    if (exists("legend_text") == FALSE) {
      legend_text <- paste0(eq, " | ", format(as.POSIXct(df_temp$time[1], format = "%Y/%m/%d %H:%M:%S"), format = "%Y/%m/%d"))
      text_cols <- cols[i]
    } else{
      legend_text <- c(legend_text, paste0(eq, " | ", format(as.POSIXct(df_temp$time[1], format = "%Y/%m/%d %H:%M:%S"), format = "%Y/%m/%d")))
      text_cols <- c(text_cols, cols[i])
    }
  }
}
legend(0, 7.25, legend = legend_text, text.col = text_cols, text.font = 2)
remove(legend_text)

#### Multi raa ####
cols <- plasma(length(input_dir) * 3)
cols <- cols[c(TRUE, FALSE, FALSE)] # Remove every other color
plot(NA, xlim = c(0, 180), ylim = c(-2.75, 6), xlab = "raa", ylab = "sif740")
for (i in 1:length(input_dir)) {
  df_temp <- get(paste0("df", i))
  if (nrow(df_temp > 0)) { # Only execute if the df has data
    points(abs(df_temp$raa), df_temp$sif740, col = cols[i])
    # polynomial line
    model <- lm(df_temp$sif740 ~ abs(df_temp$raa) + I(abs(df_temp$raa) ^ 2))
    myPredict <- predict(model)
    ix <- sort(abs(df_temp$raa), index.return = TRUE)$ix
    lines(abs(df_temp$raa)[ix], myPredict[ix], col = cols[i], lwd = 4)
    ## rounded coefficients for better output
    cf <- round(coef(model), 4)
    eq <- paste0("y = ", round(cf[1], 2),
                ifelse(sign(cf[2]) == 1, " + ", " - "), abs(round(cf[2], 2)), "x ",
                ifelse(sign(cf[3]) == 1, " + ", " - "), format(abs(cf[3]), scientific = FALSE), "x^2")
    if (exists("legend_text") == FALSE) {
      legend_text <- paste0(eq, " | ", format(as.POSIXct(df_temp$time[1], format = "%Y/%m/%d %H:%M:%S"), format = "%Y/%m/%d"))
      text_cols <- cols[i]
    } else{
      legend_text <- c(legend_text, paste0(eq, " | ", format(as.POSIXct(df_temp$time[1], format = "%Y/%m/%d %H:%M:%S"), format = "%Y/%m/%d")))
      text_cols <- c(text_cols, cols[i])
    }
  }
}
legend(0, 7.25, legend = legend_text, text.col = text_cols, text.font = 2)
remove(legend_text)

#### Single vza ####
plot(df$vza, df$sif740, xlim = c(0, 60), ylim = c(-2.75, 6), xlab = "vza", ylab = "sif740")

# polynomial line
model <- lm(df$sif740 ~ df$vza + I(df$vza ^ 2))
myPredict <- predict(model)
ix <- sort(df$vza, index.return = TRUE)$ix
lines(df$vza[ix], myPredict[ix], lwd = 4)

## rounded coefficients for better output
cf <- round(coef(model), 4)
eq <- paste0("y = ", round(cf[1], 2),
             ifelse(sign(cf[2]) == 1, " + ", " - "), abs(round(cf[2], 2)), "x ",
             ifelse(sign(cf[3]) == 1, " + ", " - "), format(abs(cf[3]), scientific = FALSE), "x^2")

legend(0, 6.75, legend = eq)

#### Single pa ####
plot(abs(df$pa), df$sif740, xlim = c(0, 70), ylim = c(-2.75, 6), xlab = "PA", ylab = "sif740")

# polynomial line
model <- lm(df$sif740 ~ abs(df$pa) + I(abs(df$pa ^ 2)))
myPredict <- predict(model)
ix <- sort(abs(df$pa), index.return = TRUE)$ix
lines(abs(df$pa[ix]), myPredict[ix], lwd = 4)

## rounded coefficients for better output
cf <- round(coef(model), 4)
eq <- paste0("y = ", round(cf[1], 2),
             ifelse(sign(cf[2]) == 1, " + ", " - "), abs(round(cf[2], 2)), "x ",
             ifelse(sign(cf[3]) == 1, " + ", " - "), format(abs(cf[3]), scientific = FALSE), "x^2")

legend(0, 6.75, legend = eq)



