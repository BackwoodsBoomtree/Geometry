files <- list.files("C:/Russell/Projects/Geometry/Data/chl", full.names = TRUE)
for(i in 1:length(files)){
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

cab_raster <- readBin("C:/Russell/Projects/Geometry/Data/chl/hd/Chlmap_phys_glb_20110625_halfDegree.bin", what = 'integer', signed = FALSE, size = 2, n = 3240000, endian = "little")
cab_raster <- matrix(data = cab_raster, nrow=360, ncol=720, byrow =T)
cab_raster <- raster(cab_raster, xmn=-180, xmx=180, ymn=-90, ymx=90, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
cab_raster[cab_raster > 150] <- NA 
plot(cab_raster)