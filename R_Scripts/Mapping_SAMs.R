library(raster)
library(rgdal)
library(ncdf4)
library(viridis)
library(leaflet)
library(mapview)

options(scipen = 999)

#### Import Data ####

file_list <- list.files(path = "G:/Geometry/oco3", full.names = TRUE, pattern = "*.nc4")

data <- nc_open(file_list[1]) # Open file

# Metadata
mode <- ncvar_get(data, "Metadata/MeasurementMode") # 0=Nadir, 1=Glint, 2=Target, 3=AreaMap, 4=Transition
orbit <- ncvar_get(data, "Metadata/OrbitId")
ID <- ncvar_get(data, "Metadata/SoundingId") # "YYYYMMDDHHMMSS"
time <- ncvar_get(data, "Delta_Time") # seconds since 1990-01-01 00:00:00 UTC
igbp <- ncvar_get(data, "Science/IGBP_index")

# SIF and Flag Data
cloud_flag <- ncvar_get(data, "Cloud/cloud_flag_abp") # 0 - \"Classified clear\", 1 - \"Classified cloudy\", 2 - \"Not classified\", all other values undefined; not used in SIF processing
q_flag <- ncvar_get(data, "Quality_Flag") # 0 = best (passes quality control + cloud fraction = 0.0); 1 = good (passes quality control); 2 = bad (failed quality control); -1 = not investigated
SIF740_D <- ncvar_get(data, "Daily_SIF_740nm")
SIF757_D <- ncvar_get(data, "Daily_SIF_757nm")
SIF771_D <- ncvar_get(data, "Daily_SIF_771nm")
SIF740 <- ncvar_get(data, "SIF_740nm")
SIF740_U <- ncvar_get(data, "SIF_Uncertainty_740nm")

# Geo Data
#lon_corners <- ncvar_get(data, "Longitude_Corners")
#lat_corners <- ncvar_get(data, "Latitude_Corners")
lon_corners <- ncvar_get(data, "Geolocation/footprint_longitude_vertices")
lat_corners <- ncvar_get(data, "Geolocation/footprint_latitude_vertices")
sza <- ncvar_get(data, "SZA")
saa <- ncvar_get(data, "SAz")
vza <- ncvar_get(data, "VZA")
vaa <- ncvar_get(data, "VAz")

# Meteo
humidity <- ncvar_get(data, "Meteo/specific_humidity")
surface_pressure <- ncvar_get(data, "Meteo/surface_pressure")
temp_skin <- ncvar_get(data, "Meteo/temperature_skin")
temp_2m <- ncvar_get(data, "Meteo/temperature_two_meter")
vpd <- ncvar_get(data, "Meteo/vapor_pressure_deficit")
wind <- ncvar_get(data, "Meteo/wind_speed")


#### Build Dataframe ####

# First, get lat/lon corners into arrays
lat1 <- lat_corners[1,]
lat2 <- lat_corners[2,]
lat3 <- lat_corners[3,]
lat4 <- lat_corners[4,]
lon1 <- lon_corners[1,]
lon2 <- lon_corners[2,]
lon3 <- lon_corners[3,]
lon4 <- lon_corners[4,]

df <- data.frame("ID" = ID, "Mode" = mode, "Cloud_Flag" = cloud_flag, "Quality_Flag" = q_flag,
                 "Time" = as.POSIXct(time, origin = "1990-01-01", tz = "UTC"),
                 "SIF740_D" = SIF740_D,"SIF757_D" = SIF757_D, "SIF771_D" = SIF771_D,
                 "SIF740" = SIF740, "SIF740_U" = SIF740_U,
                 "Lon_1" = lon1, "Lon_2" = lon2, "Lon_3" = lon3, "Lon_4" = lon4,
                 "Lat_1" = lat1, "Lat_2" = lat4, "Lat_3" = lat3, "Lat_4" = lat4,
                 "SZA" = sza, "SAA" = saa, "VZA" = vza, "VAA" = vaa,
                 "Humidity" = humidity, "Surface_Pressure" = surface_pressure,"Temp_Skin" = temp_skin,
                 "Temp_2m" = temp_2m, "VPD" = vpd, "Wind" = wind, "IGBP_Class" = igbp)

# Subset Dataframe
df_sam <- subset(df, Mode == 3)
#df_sam <- subset(df_sam, Quality_Flag == 0)
#df_sam <- subset(df_sam, Cloud_Flag == 0)

df_ATTO <- subset(df_sam, Lat_1 > 1 & Lat_1 < 4 & Lon_1 > -60 & Lon_1 < -58)

#### Create Polygons ####
for (i in 1:nrow(df_ATTO)){
  x <- c(df_ATTO$Lon_1[i],df_ATTO$Lon_2[i],df_ATTO$Lon_3[i],df_ATTO$Lon_4[i],df_ATTO$Lon_1[i])
  y <- c(df_ATTO$Lat_1[i],df_ATTO$Lat_2[i],df_ATTO$Lat_3[i],df_ATTO$Lat_4[i],df_ATTO$Lat_1[i])
  
  polly <- Polygon(cbind(x, y))
  polly <- Polygons(list(polly), row.names(df_ATTO[i, ]))
  if (i == 1){
    pollyLayer <- SpatialPolygons(list(polly))
  } else {
    pollyLayer <- SpatialPolygons(c(slot(pollyLayer, "polygons"), list(polly)))
  }
}

# Assign Data to polygons
poly_df <- SpatialPolygonsDataFrame(pollyLayer, df_ATTO)

# CRS
proj4string(poly_df) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#### Plot ####
# Get min and max lat and lon
x_min <- round(min(c(df_ATTO$Lon_1,df_ATTO$Lon_2,df_ATTO$Lon_3,df_ATTO$Lon_4)), digits=2)-0.05
x_max <- round(max(c(df_ATTO$Lon_1,df_ATTO$Lon_2,df_ATTO$Lon_3,df_ATTO$Lon_4)), digits=2)+0.05
y_min <- round(min(c(df_ATTO$Lat_1,df_ATTO$Lat_2,df_ATTO$Lat_3,df_ATTO$Lat_4)), digits=2)-0.05
y_max <- round(max(c(df_ATTO$Lat_1,df_ATTO$Lat_2,df_ATTO$Lat_3,df_ATTO$Lat_4)), digits=2)+0.05

#vPal <- viridis(nrow(df_ATTO))

### 740 ###
#df_ATTO <- df_ATTO[order(df_ATTO$SIF740),] # Sort
#df_ATTO$vPal <- vPal

#par(mfrow=c(5,5),oma=c(3,4,0.5,1))
# par(mar=c(1,1,1,1))
# par(oma=c(1,1,0.5,1))
# 
# plot(0, type='n', ann=F, axes=F, xlim=c(x_min,x_max), ylim=c(y_min,y_max))
# 
# for (i in 1:nrow(df_ATTO)){
#   polygon(c(df_ATTO$Lon_1[i],df_ATTO$Lon_2[i],df_ATTO$Lon_3[i],df_ATTO$Lon_4[i]),
#           c(df_ATTO$Lat_1[i],df_ATTO$Lat_2[i],df_ATTO$Lat_3[i],df_ATTO$Lat_4[i]),
#           col=df_ATTO$vPal[i],border=df_ATTO$vPal[i])
# }
# 
# axis(side=1, tck=0.04, mgp=c(3, 0.2, 0))
# axis(side=2, las=1, tck=0.04, mgp=c(3, 0.2, 0))
# 
# box()

vPal_rev <- colorBin(
  palette = plasma(nrow(df_ATTO)),
  domain = poly_df$SIF740,
  reverse = TRUE,
  pretty=T)

vPal <- colorBin(
  palette = plasma(nrow(df_ATTO)),
  domain = poly_df$SIF740,
  pretty=T)

par(mfrow=c(5,5),oma=c(3,4,0.5,1))
par(mar=c(1,1,1,1))
par(oma=c(1,1,0.5,1))

plot(0, type='n', ann=F, axes=F, xlim=c(x_min,x_max), ylim=c(y_min,y_max))

par(new=TRUE)

leafletOp

width = 300, height = 500,
setMaxBounds(x_min, y_min, x_max, y_max)
fitBounds(x_min, y_min, x_max, y_max)
fitBounds(x_max, y_min, x_min, y_max) %>%


# Basemap options can be found at: http://leaflet-extras.github.io/leaflet-providers/preview/index.html
m <- leaflet(poly_df, width = 550, height = 700, options = leafletOptions(zoomControl = F, attributionControl=F)) %>% 
  addProviderTiles(providers$Esri.WorldImagery) %>%
  addPolygons(data=poly_df, stroke = FALSE, smoothFactor = 0.2, fillOpacity = 1, color = ~vPal(SIF740),
              popup = paste("Sounding ID: ", poly_df$ID, "<br>",
                            "Time: ", poly_df$Time, "<br>",
                            "SIF740nm Daily: ", round(poly_df$SIF740_D,digits=2), "W/m^2/sr/�m<br>",
                            "SIF757nm Daily: ", round(poly_df$SIF757_D,digits=2), "W/m^2/sr/�m<br>",
                            "SIF771nm Daily: ", round(poly_df$SIF771_D,digits=2), "W/m^2/sr/�m<br>",
                            "SIF740nm: ", round(poly_df$SIF740,digits=2), "W/m^2/sr/�m<br>",
                            "SIF740nm Uncertainty: ", round(poly_df$SIF740_U,digits=2), "W/m^2/sr/�m<br>",
                            "Viewing Zenith Angle: ", round(poly_df$VZA,digits=2), "<br>",
                            "Viewing Azmuth Angle: ", round(poly_df$VAA,digits=2), "<br>",
                            "Solar Zenith Angle: ", round(poly_df$SZA,digits=2), "<br>",
                            "Solar Azmuth Angle: ", round(poly_df$SAA,digits=2), "<br>",
                            "Humidity: ", round(poly_df$Humidity,digits=4), "kg/kg<br>",
                            "Surface Pressure: ", round(poly_df$Surface_Pressure,digits=2), "Pa<br>",
                            "Temperature Skin: ", round(poly_df$Temp_Skin,digits=2), "K<br>",
                            "Temperature 2m: ", round(poly_df$Temp_2m,digits=2), "K<br>",
                            "VPD: ", round(poly_df$VPD,digits=2), "Pa<br>",
                            "Wind Speed: ", round(poly_df$Wind,digits=2), "m/s<br>",
                            "IGBP Class: ", poly_df$IGBP_Class))  %>%
  addLegend("bottomright", pal = vPal_rev, values = ~SIF740, title = "SIF740nm<br>W/m^2/sr/�m", opacity = 1,
            labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))) %>%
  addMiniMap("topright", tiles = providers$CartoDB.Voyager, width = 100, height = 100,
             toggleDisplay = F, zoomLevelFixed = 1,
             aimingRectOptions = list(color = "#d91414", weight  = 3, clickable = FALSE),
             shadowRectOptions = list(color = "#000000", weight = 1, clickable = FALSE, opacity = 0, fillOpacity = 0))

m

mapshot(m, file = "C:/Russell/R_Scripts/Geometry/Mapping_SAMs.pdf")

par(mfrow=c(5,5),oma=c(3,4,0.5,1))
par(mar=c(1,1,1,1))
par(oma=c(1,1,0.5,1))

plot(m, type='n', ann=F, axes=F, xlim=c(x_min,x_max), ylim=c(y_min,y_max))

