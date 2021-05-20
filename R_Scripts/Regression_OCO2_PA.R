library(dplyr)
library(tidyr)
library(lubridate)
library(LSD)
library(MASS)

# PA  <- abs(df111$PA)
# SIF <- df111$SIF_757nm

PA  <- c(1, 2, 3, 4, 5, 6)
SIF <- c(-1, 1, 0, 1.25, 1.5, 2)

par(oma=c(3,3,2,1.25))

############# Plot #####################

# Axes margin
op <- par(mar = c(0,0,0,0))

heatscatter(PA, SIF, cex.axis=1.25, ylim=c(-1,2), xlim=c(0,70),
            cexplot = 1, main=FALSE, tck=0.03, mgp=c(3, 0.2, 0))

# polynomial line
model <- lm(SIF ~ PA+I(PA^2))
myPredict <- predict(model)
ix <- sort(PA,index.return=T)$ix
lines(PA[ix], myPredict[ix], col=rgb(127,188,65,max=255), lwd=4)

#text(35,45, "C", cex=1.25)

## rounded coefficients for better output
cf <- round(coef(model), 4) 

## sign check to avoid having plus followed by minus for negative coefficients
eq <- paste0("y = ", round(cf[1], 2),
             ifelse(sign(cf[2])==1, " + ", " - "), abs(round(cf[2],2)), "x ",
             ifelse(sign(cf[3])==1, " + ", " - "), format(abs(cf[3]), scientific=FALSE), "x^2")

## printing of the equation
title(eq, line=-1)

box()

############################

mtext(3, text = "Lamont, Oklahoma 2015-05-17", cex=1.25)
mtext(2, text=expression(paste("OCO-2 SIF"['757']*" (mW/m"^{2}*"/sr/nm)")), line = 1.25, cex=1.25)
mtext(1, text=expression(paste("Phase Angle (absolute)")), line=2, cex=1.25)
