library(raster)
library(rgdal)
library(viridis)

#### Load Map ####

coastlines <- readOGR("C:/Russell/R_Scripts/TROPOMI_2/mapping/GSHHS_shp/c/GSHHS_c_L1.shp")
class(coastlines)
extent(coastlines)
crs(coastlines)

# Data
raster_gosat_annual_2019 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/GOSAT_Annual_Mean_757Daily_2019_0.50.tif")
raster_gosat_june_2019   <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/GOSAT_June_Mean_757Daily_2019_0.50.tif")
raster_oco2_annual_2020  <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/OCO2_Annual_Mean_757Daily_2020_0.50.tif")
raster_oco2_june_2020    <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/OCO2_June_Mean_757Daily_2020_0.50.tif")
raster_oco3_annual_2020  <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/OCO3_Annual_Mean_757Daily_2020_0.50.tif")
raster_oco3_june_2020    <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/OCO3_June_Mean_757Daily_2020_0.50.tif")

# Rescale for visualization
raster_gosat_annual_2019[raster_gosat_annual_2019 < 0]   <- 0
raster_gosat_annual_2019[raster_gosat_annual_2019 > 0.8] <- 0.8
raster_oco2_annual_2020[raster_oco2_annual_2020 < 0]     <- 0
raster_oco2_annual_2020[raster_oco2_annual_2020 > 0.8]   <- 0.8
raster_oco3_annual_2020[raster_oco3_annual_2020 < 0]     <- 0
raster_oco3_annual_2020[raster_oco3_annual_2020 > 0.8]   <- 0.8

raster_gosat_june_2019[raster_gosat_june_2019 < 0]   <- 0
raster_gosat_june_2019[raster_gosat_june_2019 > 0.8] <- 0.8
raster_oco2_june_2020[raster_oco2_june_2020 < 0]     <- 0
raster_oco2_june_2020[raster_oco2_june_2020 > 0.8]   <- 0.8
raster_oco3_june_2020[raster_oco3_june_2020 < 0]     <- 0
raster_oco3_june_2020[raster_oco3_june_2020 > 0.8]   <- 0.8

# Colors
sif.col <- viridis(8)

pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Grid_GOSAT_OCO/Plot_Annual_Monthly_GOSAT_OCO.pdf", width = 7.5, height = 5.75,compress=FALSE)

par(mfrow=c(3,2),oma=c(1.5,0.25,1.75,0))

## GOSAT Annual ##
op <- par(mar = c(0,0,0,0.25))
plot(raster_gosat_annual_2019, col = sif.col,  ext=c(-180, 180, -60, 85), axes = F, xaxs = "i", yaxs = "i", horizontal = T, legend = F)
mtext(3, text = expression(paste("GOSAT 2019 Annual Mean Daily SIF"['757']*"")))
mtext(3, text = "A", font = 2, at = -170, line = -1.5)
plot(coastlines, add = TRUE, lwd = 0.5)

# plot(Nmean, legend.only=TRUE, col=N.cols, horizontal=F, legend.width=2, legend.shrink=0.75,
#      legend.args = list(text=expression(paste("Retrievals")), side = 3, line = 0.25, cex=0.85),
#      axis.args = list(line = -0.75, cex.axis=0.85,tick=F, at=c(1,66), labels=c("1","66")),
#      smallplot=c(0.075,0.12,0.37,0.6)); par(mar = par("mar"))

## GOSAT June ##
op <- par(mar = c(0,0,0,0.25))
plot(raster_gosat_june_2019, col = sif.col,  ext=c(-180, 180, -60, 85), axes = F, xaxs = "i", yaxs = "i", horizontal = T, legend = F)
mtext(3, text = expression(paste("GOSAT June 2019 Mean Daily SIF"['757']*"")))
mtext(3, text = "B", font = 2, at = -170, line = -1.5)
plot(coastlines, add = TRUE, lwd = 0.5)

## OCO2 Annual ##
op <- par(mar = c(0,0,0,0.25))
plot(raster_oco2_annual_2020, col = sif.col,  ext=c(-180, 180, -60, 85), axes = F, xaxs = "i", yaxs = "i", horizontal = T, legend = F)
mtext(3, text = expression(paste("OCO-2 2020 Annual Mean Daily SIF"['757']*"")))
mtext(3, text = "C", font = 2, at = -170, line = -1.5)
plot(coastlines, add = TRUE, lwd = 0.5)

## OCO2 June ##
op <- par(mar = c(0,0,0,0.25))
plot(raster_oco2_june_2020, col = sif.col,  ext=c(-180, 180, -60, 85), axes = F, xaxs = "i", yaxs = "i", horizontal = T, legend = F)
mtext(3, text = expression(paste("OCO-2 June 2020 Mean Daily SIF"['757']*"")))
mtext(3, text = "D", font = 2, at = -170, line = -1.5)
plot(coastlines, add = TRUE, lwd = 0.5)

## OCO3 Annual ##
op <- par(mar = c(0,0,0,0.25))
plot(raster_oco3_annual_2020, col = sif.col,  ext=c(-180, 180, -60, 85), axes = F, xaxs = "i", yaxs = "i", horizontal = T, legend = F)
mtext(3, text = expression(paste("OCO-3 2020 Annual Mean Daily SIF"['757']*"")))
mtext(3, text = "E", font = 2, at = -170, line = -1.5)
plot(coastlines, add = TRUE, lwd = 0.5)

## OCO3 June ##
op <- par(mar = c(0,0,0,0.25))
plot(raster_oco3_june_2020, col = sif.col,  ext=c(-180, 180, -60, 85), axes = F, xaxs = "i", yaxs = "i", horizontal = T, legend = F)
mtext(3, text = expression(paste("OCO-3 June 2020 Mean Daily SIF"['757']*"")))
mtext(3, text = "F", font = 2, at = -170, line = -1.5)
plot(coastlines, add = TRUE, lwd = 0.5)

# plot(raster_oco3_june_2020, legend.only = TRUE, col = sif.col, horizontal = T, legend.width = 2, legend.shrink = 0.75,
#      legend.args = list(text=expression(paste("W/m"^"-2"*"/sr/µm")), side = 1, line = 0.2, cex = 0.85),
#      axis.args = list(line = -1.05, cex.axis = 1, tick = FALSE, at = c(0, 0.8), labels = c("<0", ">0.8")),
#      smallplot=c(0.40, 0.90, 0.265, 0.315)); par(mar = par("mar"))

# Move in Adobe PDF
plot(raster_oco3_june_2020, legend.only = TRUE, col = sif.col, horizontal = T, legend.width = 2, legend.shrink = 0.75,
     legend.args = list(text=expression(paste("W/m"^"-2"*"/sr/µm")), side = 1, line = 1, cex = 1),
     axis.args = list(line = -0.50, cex.axis = 1.25, tick = FALSE, at = c(0, 0.8), labels = c("<0", ">0.8")))


dev.off()
