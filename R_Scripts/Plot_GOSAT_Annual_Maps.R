library(raster)
library(rgdal)
library(RColorBrewer)
library(viridis)

raster_gosat_annual_2015 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2015.tif")
raster_gosat_annual_2016 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2016.tif")
raster_gosat_annual_2017 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2017.tif")
raster_gosat_annual_2018 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2018.tif")
raster_gosat_annual_2019 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2019.tif")

# Set thresholds for visual mapping
raster_gosat_annual_t_2015 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2015.tif")
raster_gosat_annual_t_2016 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2016.tif")
raster_gosat_annual_t_2017 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2017.tif")
raster_gosat_annual_t_2018 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2018.tif")
raster_gosat_annual_t_2019 <- raster("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean_740Daily_2019.tif")

raster_gosat_annual_t_2015[raster_gosat_annual_t_2015 < 0] <- 0
raster_gosat_annual_t_2015[raster_gosat_annual_t_2015 > 0.7]  <- 0.7
raster_gosat_annual_t_2016[raster_gosat_annual_t_2016 < 0] <- 0
raster_gosat_annual_t_2016[raster_gosat_annual_t_2016 > 0.7]  <- 0.7
raster_gosat_annual_t_2017[raster_gosat_annual_t_2017 < 0] <- 0
raster_gosat_annual_t_2017[raster_gosat_annual_t_2017 > 0.7]  <- 0.7
raster_gosat_annual_t_2018[raster_gosat_annual_t_2018 < 0] <- 0
raster_gosat_annual_t_2018[raster_gosat_annual_t_2018 > 0.7]  <- 0.7
raster_gosat_annual_t_2019[raster_gosat_annual_t_2019 < 0] <- 0
raster_gosat_annual_t_2019[raster_gosat_annual_t_2019 > 0.7]  <- 0.7

round2 <- function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

labs <- c(expression(paste("Mean GOSAT Daily SIF"['740']*" 2015")),
          expression(paste("Mean GOSAT Daily SIF"['740']*" 2016")),
          expression(paste("Mean GOSAT Daily SIF"['740']*" 2017")),
          expression(paste("Mean GOSAT Daily SIF"['740']*" 2018")),
          expression(paste("Mean GOSAT Daily SIF"['740']*" 2019")))

### Plotting ###

coastlines <- readOGR("C:/Russell/R_Scripts/TROPOMI_2/mapping/GSHHS_shp/c/GSHHS_c_L1.shp")
class(coastlines)
extent(coastlines)
crs(coastlines)

sif.col <- viridis(8)

pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/GOSAT_Annual_Mean.pdf", width=7, height=4.25, compress=FALSE)

par(mfrow=c(3, 3), oma=c(0, 0.25, 1.25, 0))

### GOSAT 2015 ###
meanSIF <- round2(cellStats(raster_gosat_annual_2015, mean, na.rm=T), 2)
op <- par(mar = c(0, 0, 0.25, 0.25))
plot(raster_gosat_annual_t_2015, col = sif.col, ext = c(-180, 180, -56, 80), axes = FALSE, xaxs = "i", yaxs = "i", horizontal = TRUE, legend = FALSE)
mtext(3, text=labs[1], cex = 0.85)
mtext(3, text = "A", cex = 0.85, adj = 0, font = 2)
plot(coastlines, add = TRUE, lwd=0.5)

plot(raster_gosat_annual_t_2015, legend.only = TRUE, col = sif.col, horizontal = T, legend.width = 2, legend.shrink = 0.75,
     legend.args = list(text=expression(paste("W/m"^"-2"*"/sr/µm")), side = 1, line = 0.2, cex = 0.60),
     axis.args = list(line = -1.05, cex.axis = 1, tick = FALSE, at = c(0, 0.7), labels = c("<0", ">0.7")),
     smallplot=c(0.40, 0.90, 0.265, 0.315)); par(mar = par("mar"))

par(new = TRUE)
op <- par(mar = c(3.1, 0.25, 5, 13)) # Set margins
hist(raster_gosat_annual_2015, col = "gray75", breaks = 30, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white", border = NA)
abline(v = meanSIF, col = "red")
abline(h = 100000, lty = 2)
axis(3, tck = FALSE, labels = meanSIF, at = meanSIF, mgp = c(3, 0.1, 0))
axis(1, tck = FALSE, labels = c(-1, 0, 1), at = c(-1, 0, 1), mgp = c(3, 0.1, 0))
hist(raster_gosat_annual_2015, col = "gray75", breaks = 30, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE, add = T)
box()

### GOSAT 2016 ###
meanSIF <- round2(cellStats(raster_gosat_annual_2016, mean, na.rm=T), 2)
op <- par(mar = c(0, 0, 0.25, 0.25))
plot(raster_gosat_annual_t_2016, col = sif.col, ext = c(-180, 180, -56, 80), axes = FALSE, xaxs = "i", yaxs = "i", horizontal = TRUE, legend = FALSE)
mtext(3, text=labs[2], cex = 0.85)
mtext(3, text = "B", cex = 0.85, adj = 0, font = 2)
plot(coastlines, add = TRUE, lwd=0.5)

plot(raster_gosat_annual_t_2016, legend.only = TRUE, col = sif.col, horizontal = T, legend.width = 2, legend.shrink = 0.75,
     legend.args = list(text=expression(paste("W/m"^"-2"*"/sr/µm")), side = 1, line = 0.2, cex = 0.60),
     axis.args = list(line = -1.05, cex.axis = 1, tick = FALSE, at = c(0, 0.7), labels = c("<0", ">0.7")),
     smallplot=c(0.40, 0.90, 0.265, 0.315)); par(mar = par("mar"))

par(new = TRUE)
op <- par(mar = c(3.1, 0.25, 5, 13)) # Set margins
hist(raster_gosat_annual_2016, col = "gray75", breaks = 30, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white", border = NA)
abline(v = meanSIF, col = "red")
abline(h = 100000, lty = 2)
axis(3, tck = FALSE, labels = meanSIF, at = meanSIF, mgp = c(3, 0.1, 0))
axis(1, tck = FALSE, labels = c(-1, 0, 1), at = c(-1, 0, 1), mgp = c(3, 0.1, 0))
hist(raster_gosat_annual_2016, col = "gray75", breaks = 30, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE, add = T)
box()

### GOSAT 2017 ###
meanSIF <- round2(cellStats(raster_gosat_annual_2017, mean, na.rm=T), 2)
op <- par(mar = c(0, 0, 0.25, 0.25))
plot(raster_gosat_annual_t_2017, col = sif.col, ext = c(-180, 180, -56, 80), axes = FALSE, xaxs = "i", yaxs = "i", horizontal = TRUE, legend = FALSE)
mtext(3, text=labs[3], cex = 0.85)
mtext(3, text = "C", cex = 0.85, adj = 0, font = 2)
plot(coastlines, add = TRUE, lwd=0.5)

plot(raster_gosat_annual_t_2017, legend.only = TRUE, col = sif.col, horizontal = T, legend.width = 2, legend.shrink = 0.75,
     legend.args = list(text=expression(paste("W/m"^"-2"*"/sr/µm")), side = 1, line = 0.2, cex = 0.60),
     axis.args = list(line = -1.05, cex.axis = 1, tick = FALSE, at = c(0, 0.7), labels = c("<0", ">0.7")),
     smallplot=c(0.40, 0.90, 0.265, 0.315)); par(mar = par("mar"))

par(new = TRUE)
op <- par(mar = c(3.1, 0.25, 5, 13)) # Set margins
hist(raster_gosat_annual_2017, col = "gray75", breaks = 30, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white", border = NA)
abline(v = meanSIF, col = "red")
abline(h = 100000, lty = 2)
axis(3, tck = FALSE, labels = meanSIF, at = meanSIF, mgp = c(3, 0.1, 0))
axis(1, tck = FALSE, labels = c(-1, 0, 1), at = c(-1, 0, 1), mgp = c(3, 0.1, 0))
hist(raster_gosat_annual_2017, col = "gray75", breaks = 30, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE, add = T)
box()

### GOSAT 2018 ###
meanSIF <- round2(cellStats(raster_gosat_annual_2018, mean, na.rm=T), 2)
op <- par(mar = c(0, 0, 0.25, 0.25))
plot(raster_gosat_annual_t_2018, col = sif.col, ext = c(-180, 180, -56, 80), axes = FALSE, xaxs = "i", yaxs = "i", horizontal = TRUE, legend = FALSE)
mtext(3, text=labs[4], cex = 0.85)
mtext(3, text = "D", cex = 0.85, adj = 0, font = 2)
plot(coastlines, add = TRUE, lwd=0.5)

plot(raster_gosat_annual_t_2018, legend.only = TRUE, col = sif.col, horizontal = T, legend.width = 2, legend.shrink = 0.75,
     legend.args = list(text=expression(paste("W/m"^"-2"*"/sr/µm")), side = 1, line = 0.2, cex = 0.60),
     axis.args = list(line = -1.05, cex.axis = 1, tick = FALSE, at = c(0, 0.7), labels = c("<0", ">0.7")),
     smallplot=c(0.40, 0.90, 0.265, 0.315)); par(mar = par("mar"))

par(new = TRUE)
op <- par(mar = c(3.1, 0.25, 5, 13)) # Set margins
hist(raster_gosat_annual_2018, col = "gray75", breaks = 20, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white", border = NA)
abline(v = meanSIF, col = "red")
abline(h = 100000, lty = 2)
axis(3, tck = FALSE, labels = meanSIF, at = meanSIF, mgp = c(3, 0.1, 0))
axis(1, tck = FALSE, labels = c(-1, 0, 1), at = c(-1, 0, 1), mgp = c(3, 0.1, 0))
hist(raster_gosat_annual_2018, col = "gray75", breaks = 20, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE, add = T)
box()

### GOSAT 2019 ###
meanSIF <- round2(cellStats(raster_gosat_annual_2019, mean, na.rm=T), 2)
op <- par(mar = c(0, 0, 0.25, 0.25))
plot(raster_gosat_annual_t_2019, col = sif.col, ext = c(-180, 180, -56, 80), axes = FALSE, xaxs = "i", yaxs = "i", horizontal = TRUE, legend = FALSE)
mtext(3, text=labs[5], cex = 0.85)
mtext(3, text = "E", cex = 0.85, adj = 0, font = 2)
plot(coastlines, add = TRUE, lwd=0.5)

plot(raster_gosat_annual_t_2019, legend.only = TRUE, col = sif.col, horizontal = T, legend.width = 2, legend.shrink = 0.75,
     legend.args = list(text=expression(paste("W/m"^"-2"*"/sr/µm")), side = 1, line = 0.2, cex = 0.60),
     axis.args = list(line = -1.05, cex.axis = 1, tick = FALSE, at = c(0, 0.7), labels = c("<0", ">0.7")),
     smallplot=c(0.40, 0.90, 0.265, 0.315)); par(mar = par("mar"))

par(new = TRUE)
op <- par(mar = c(3.1, 0.25, 5, 13)) # Set margins
hist(raster_gosat_annual_2019, col = "gray75", breaks = 30, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white", border = NA)
abline(v = meanSIF, col = "red")
abline(h = 100000, lty = 2)
axis(3, tck = FALSE, labels = meanSIF, at = meanSIF, mgp = c(3, 0.1, 0))
axis(1, tck = FALSE, labels = c(-1, 0, 1), at = c(-1, 0, 1), mgp = c(3, 0.1, 0))
hist(raster_gosat_annual_2019, col = "gray75", breaks = 30, ylim = c(0, 1000), xlim = c(-1, 1), xaxs = "i", yaxs = "i", ann = FALSE, axes = FALSE, add = T)
box()

dev.off()