library(raster)
library(ncdf4)

# Create Weekly Means
folders <- list.dirs("C:/Users/rusty/Downloads/Leaf_Chlorophyll_Map_0.5deg", full.names = TRUE, recursive = FALSE)
for (i in 1:52){
  for (j in 1:length(folders)){
    cab_raster <- list.files(folders[j], full.names = TRUE)[i]
    cab_raster <- readBin(cab_raster, what = 'integer', signed = FALSE, size = 2, n = 3240000, endian = "little")
    cab_raster <- matrix(data = cab_raster, nrow=360, ncol=720, byrow =T)
    cab_raster <- raster(cab_raster, xmn=-180, xmx=180, ymn=-90, ymx=90, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    cab_raster[cab_raster > 150] <- NA
    if (j == 1){
      cab_brick <- brick(cab_raster)
    } else {
      cab_brick <- addLayer(cab_brick, cab_raster)
    }
  }
  cab_mean <- calc(cab_brick, fun = mean, na.rm = TRUE)
  writeRaster(cab_mean, filename = paste0("C:/Users/rusty/Downloads/Cab_mean/Chl_Mean_2003-2011_Week_", sprintf("%02d", i), ".tif"), overwrite = TRUE)
}

# Create NC file of weekly means
files <- list.files("C:/Users/rusty/Downloads/Cab_mean", full.names = TRUE)
for (i in 1:length(files)){
  if (i == 1){
    cab_brick <- brick(files[i])
  } else {
    cab_brick <- addLayer(cab_brick, files[i])
  }
}

writeRaster(cab_brick, "C:/Users/rusty/Downloads/Cab_mean/Chl_Mean_2003-2011_Weekly_0.05deg.nc", overwrite = TRUE, format = "CDF", varname = "chl", varunit = "μg cm-2",
        longname = "Leaf Chlorophyll Content", xname = "lon", yname = "lat", zname = "time", zunit = "Week", options="COMPRESS=LZW")




### Exploring 300-m data
files <- list.files("C:/Russell/Projects/Geometry/Data/chl", full.names = TRUE)
for (i in 1:length(files)){
  cab_raster <- readBin(files[12], what = 'integer', signed = FALSE, size = 2, n = 3240000, endian = "little")
  cab_raster <- matrix(data = cab_raster, nrow = 1800, ncol = 1800, byrow = TRUE)
  cab_raster <- raster(cab_raster, xmn = -60, xmx = -55, ymn = 0, ymx = 5, crs = "+init=epsg:4326")
  cab_raster[cab_raster > 150] <- NA    # optional to screen out flagged data for visual assessment
  if(i == 1){
    cab_brick <- brick(cab_raster)
    } else {
      cab_brick <- addLayer(cab_brick, cab_raster)
  }
}
cab_mean <- calc(cab_brick, fun = mean, na.rm = TRUE)

cab_raster <- readBin("C:/Russell/Projects/Geometry/Data/chl/Chlmap_phys-GFwith2010-h24v17-20111224_lacc.bin", what = 'integer', signed = FALSE, size = 2, n = 3240000, endian = "little")
cab_raster <- matrix(data = cab_raster, nrow = 1800, ncol = 1800, byrow = TRUE)
cab_raster <- raster(cab_raster, xmn = -60, xmx = -55, ymn = 0, ymx = 5, crs = "+init=epsg:4326")
cab_raster[cab_raster > 150] <- NA    # optional to screen out flagged data for visual assessment
spplot(cab_raster)
plot(cab_raster)