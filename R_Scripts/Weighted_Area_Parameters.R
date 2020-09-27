library(sp)
library(ncdf4)
library(raster)
library(rgdal)

lai <- raster("C:/Russell/Projects/Geometry/Data/lai/c_gls_LAI300-LAI_201906200000_CUSTOM_PROBAV_V1.0.1.tiff")

poly_df_wgs84 <- poly_df
crs(poly_df_wgs84) <- crs(lai)

lai <- crop(lai, poly_df_wgs84)

# Sample parameter raster using OCO footprint polygons and calculate area weighted mean
lai_awm <- extract(lai, poly_df, weights = TRUE, fun = mean, na.rm = TRUE)

poly_df_wgs84$lai <- as.vector(lai_awm)

# Export df to csv
write.csv(as.data.frame(poly_df_wgs84), paste0("C:/Russell/R_Scripts/Geometry/sif_ATTO_Tower_Manaus_Brazil_(incorrect)_", format(lab_time, format = '%Y-%m-%d', usetz = FALSE), ".csv"), row.names = FALSE)

# Export to shapefile
shapefile(poly_df_wgs84, paste0("C:/Russell/R_Scripts/Geometry/sif_ATTO_Tower_Manaus_Brazil_(incorrect)_", format(lab_time, format = '%Y-%m-%d', usetz = FALSE), ".shp"))

View(as.data.frame(poly_df_wgs84))

pdf("C:/Russell/R_Scripts/Geometry/test.pdf")
test <- plot(lai)
test <- plot(poly_df, add = TRUE)
dev.off()