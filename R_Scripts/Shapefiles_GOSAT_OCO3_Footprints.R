library(raster)
library(rgdal)
library(ncdf4)


unfiltered_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/unfiltered_matched_sounding_means.csv")
veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/veg_matched_sounding_means.csv")
time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/time_veg_matched_sounding_means.csv")
vza_time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/vza_time_veg_matched_sounding_means.csv")
n_vza_time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/n_vza_time_veg_matched_sounding_means.csv")
temp_n_vza_time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/temp_n_vza_time_veg_matched_sounding_means.csv")

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

#region ############# Plotting overlaps as points ############

# All soundings unfiltered
coords   <- as.data.frame(cbind(unfiltered_matched_sounding_means$longitude, unfiltered_matched_sounding_means$latitude))
colnames(coords) <- c("longitude", "latitude")
point_df <- SpatialPointsDataFrame(coords, proj4string = CRS("+init=epsg:4326"), data = unfiltered_matched_sounding_means, coords.nrs = c(8, 9))
shapefile(point_df, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/shp/unfiltered_matched_sounding_means.shp", overwrite = TRUE)

# Soundings with all filters
coords   <- as.data.frame(cbind(temp_n_vza_time_veg_matched_sounding_means$longitude, temp_n_vza_time_veg_matched_sounding_means$latitude))
colnames(coords) <- c("longitude", "latitude")
point_df <- SpatialPointsDataFrame(coords, proj4string = CRS("+init=epsg:4326"), data = temp_n_vza_time_veg_matched_sounding_means, coords.nrs = c(8, 9))
shapefile(point_df, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/shp/temp_n_vza_time_veg_matched_sounding_means.shp", overwrite = TRUE)

# Soundings filtered only for vegetation
# Soundings with all filters
coords   <- as.data.frame(cbind(veg_matched_sounding_means$longitude, veg_matched_sounding_means$latitude))
colnames(coords) <- c("longitude", "latitude")
point_df <- SpatialPointsDataFrame(coords, proj4string = CRS("+init=epsg:4326"), data = veg_matched_sounding_means, coords.nrs = c(8, 9))
shapefile(point_df, "C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/shp/veg_matched_sounding_means.shp", overwrite = TRUE)

###############################

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