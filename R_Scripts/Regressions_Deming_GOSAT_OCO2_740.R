library(mcr)


unfiltered_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/unfiltered_matched_sounding_means.csv")
veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/veg_matched_sounding_means.csv")
time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/time_veg_matched_sounding_means.csv")
vza_time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/vza_time_veg_matched_sounding_means.csv")
n_vza_time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/n_vza_time_veg_matched_sounding_means.csv")
temp_n_vza_time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/temp_n_vza_time_veg_matched_sounding_means.csv")


#### Plot titles ####

title_1 <- "QC = Best; Cloud Free"
title_2 <- "Vegetation Only"
title_3 <- "Veg;  < 1 hour"
title_4 <- "Veg; < 1 hour; VZA < 5"
title_5 <- "Veg; < 1 hour; VZA < 5; N >= 10"
title_6 <- "Veg; < 1 hour; VZA < 5; N >= 10;\nT >= 5C"

#region ############# PLOTS WITH DIFFERENT FILTERS - GOSAT P #############

### STATS ###

# Mean up the GOSAT polarizations
mean_GOSAT_SIF_740nm_unfiltered          <- unfiltered_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_veg                 <- veg_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_time_veg            <- time_veg_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_vza_time_veg        <- vza_time_veg_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_n_vza_time_veg      <- n_vza_time_veg_matched_sounding_means$SIF_740nm_P
mean_GOSAT_SIF_740nm_temp_n_vza_time_veg <- temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_P

# Variances of the errors
var_err_GOSAT_SIF_740nm_unfiltered          <- var(unfiltered_matched_sounding_means$SIF_Uncertainty_740nm_P)
var_err_GOSAT_SIF_740nm_veg                 <- var(veg_matched_sounding_means$SIF_Uncertainty_740nm_P)
var_err_GOSAT_SIF_740nm_time_veg            <- var(time_veg_matched_sounding_means$SIF_Uncertainty_740nm_P)
var_err_GOSAT_SIF_740nm_vza_time_veg        <- var(vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_P)
var_err_GOSAT_SIF_740nm_n_vza_time_veg      <- var(n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_P)
var_err_GOSAT_SIF_740nm_temp_n_vza_time_veg <- var(temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_P)

var_err_OCO2_SIF_740nm_unfiltered          <- var(unfiltered_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_veg                 <- var(veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_time_veg            <- var(time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_vza_time_veg        <- var(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_n_vza_time_veg      <- var(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_temp_n_vza_time_veg <- var(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)

var_err_ratio_SIF_740nm_unfiltered          <- var_err_OCO2_SIF_740nm_unfiltered / var_err_GOSAT_SIF_740nm_unfiltered
var_err_ratio_SIF_740nm_veg                 <- var_err_OCO2_SIF_740nm_veg / var_err_GOSAT_SIF_740nm_veg
var_err_ratio_SIF_740nm_time_veg            <- var_err_OCO2_SIF_740nm_time_veg / var_err_GOSAT_SIF_740nm_time_veg 
var_err_ratio_SIF_740nm_vza_time_veg        <- var_err_OCO2_SIF_740nm_vza_time_veg / var_err_GOSAT_SIF_740nm_vza_time_veg
var_err_ratio_SIF_740nm_n_vza_time_veg      <- var_err_OCO2_SIF_740nm_n_vza_time_veg / var_err_GOSAT_SIF_740nm_n_vza_time_veg
var_err_ratio_SIF_740nm_temp_n_vza_time_veg <- var_err_OCO2_SIF_740nm_temp_n_vza_time_veg / var_err_GOSAT_SIF_740nm_temp_n_vza_time_veg

# Run Deming regressions
reg_SIF_740nm_unfiltered          <- mcreg(mean_GOSAT_SIF_740nm_unfiltered, unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_unfiltered,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_veg                 <- mcreg(mean_GOSAT_SIF_740nm_veg, veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_time_veg            <- mcreg(mean_GOSAT_SIF_740nm_time_veg, time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_vza_time_veg        <- mcreg(mean_GOSAT_SIF_740nm_vza_time_veg, vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_vza_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_n_vza_time_veg      <- mcreg(mean_GOSAT_SIF_740nm_n_vza_time_veg, n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_n_vza_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_temp_n_vza_time_veg <- mcreg(mean_GOSAT_SIF_740nm_temp_n_vza_time_veg, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_temp_n_vza_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")

# Run linear regressions for R2 values
lm_SIF_740nm_unfiltered          <- lm(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered)
lm_SIF_740nm_veg                 <- lm(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg)
lm_SIF_740nm_time_veg            <- lm(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg)
lm_SIF_740nm_vza_time_veg        <- lm(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg)
lm_SIF_740nm_n_vza_time_veg      <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg)
lm_SIF_740nm_temp_n_vza_time_veg <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg)

# Summaries
summary_unfiltered          <- summary(lm_SIF_740nm_unfiltered)
summary_veg                 <- summary(lm_SIF_740nm_veg)
summary_time_veg            <- summary(lm_SIF_740nm_time_veg)
summary_vza_time_veg        <- summary(lm_SIF_740nm_vza_time_veg)
summary_n_vza_time_veg      <- summary(lm_SIF_740nm_n_vza_time_veg)
summary_temp_n_vza_time_veg <- summary(lm_SIF_740nm_temp_n_vza_time_veg)

# Round coefficients
cf_lm_SIF_740nm_unfiltered           <- round(coef(lm_SIF_740nm_unfiltered), 2)
cf_lm_SIF_740nm_veg                  <- round(coef(lm_SIF_740nm_veg), 2)
cf_lm_SIF_740nm_time_veg             <- round(coef(lm_SIF_740nm_time_veg), 2)
cf_lm_SIF_740nm_vza_time_veg         <- round(coef(lm_SIF_740nm_vza_time_veg), 2)
cf_lm_SIF_740nm_n_vza_time_veg       <- round(coef(lm_SIF_740nm_n_vza_time_veg), 2)
cf_lm_SIF_740nm_temp_n_vza_time_veg  <- round(coef(lm_SIF_740nm_temp_n_vza_time_veg), 2)

# Equations
eq_unfiltered          <- paste0("y = ", round(reg_SIF_740nm_unfiltered@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_unfiltered@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_unfiltered@para[2], 2)), "x")
eq_veg                 <- paste0("y = ", round(reg_SIF_740nm_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_veg@para[2], 2)), "x")
eq_time_veg            <- paste0("y = ", round(reg_SIF_740nm_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_time_veg@para[2], 2)), "x")
eq_vza_time_veg        <- paste0("y = ", round(reg_SIF_740nm_vza_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_vza_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_vza_time_veg@para[2], 2)), "x")
eq_n_vza_time_veg      <- paste0("y = ", round(reg_SIF_740nm_n_vza_time_veg@para[1], 2),
                                ifelse(sign(round(reg_SIF_740nm_n_vza_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_n_vza_time_veg@para[2], 2)), "x")
eq_temp_n_vza_time_veg <- paste0("y = ", round(reg_SIF_740nm_temp_n_vza_time_veg@para[1], 2),
                                ifelse(sign(round(reg_SIF_740nm_temp_n_vza_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_temp_n_vza_time_veg@para[2], 2)), "x")

# R values
r_unfiltered          <- paste0("p < .001, R2 = ", round(summary_unfiltered$adj.r.squared, 2))
r_veg                 <- paste0("p < .001, R2 = ", round(summary_veg$adj.r.squared, 2))
r_time_veg            <- paste0("p < .001, R2 = ", round(summary_time_veg$adj.r.squared, 2))
r_vza_time_veg        <- paste0("p < .001, R2 = ", round(summary_vza_time_veg$adj.r.squared, 2))
r_n_vza_time_veg      <- paste0("p < .001, R2 = ", round(summary_n_vza_time_veg$adj.r.squared, 2))
r_temp_n_vza_time_veg <- paste0("p < .001, R2 = ", round(summary_temp_n_vza_time_veg$adj.r.squared, 2))


##### SETUP MAIN PLOT #####

pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_Deming_GOSAT_P_OCO2_SIF740.pdf", width=8.5, height=6, compress=FALSE)

par(mfrow=c(2, 3), oma = c(0, 1.5, 0, 3.5))

# FILTERED DATA PLOTS

# Unfiltered
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_unfiltered, x.lab = "", y.lab = "", main = title_1, xlim = c(-6, 6), ylim = c(-6, 6),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_unfiltered, r_unfiltered), bty = 'o', bg = "white", box.col = "white")

# Veg
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_veg, x.lab = "", y.lab = "", main = title_2, xlim = c(-6, 6), ylim = c(-6, 6),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_veg, r_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_time_veg, x.lab = "", y.lab = "", main = title_3, xlim = c(-6, 6), ylim = c(-6, 6),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_time_veg, r_time_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time VZA
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_vza_time_veg, x.lab = "", y.lab = "", main = title_4, xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_vza_time_veg, r_vza_time_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time VZA N
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_n_vza_time_veg, x.lab = "", y.lab = "", main = title_5, xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_n_vza_time_veg, r_n_vza_time_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time VZA N Temp
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_temp_n_vza_time_veg, x.lab = "", y.lab = "", main = title_6, xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_temp_n_vza_time_veg, r_temp_n_vza_time_veg), bty = 'o', bg = "white", box.col = "white")

mtext(expression(paste("Mean OCO SIF 740nm")), 2, -0.5, cex = 1, outer = TRUE)
mtext(expression(paste("GOSAT P SIF 740nm")), 1, -1.5, cex = 1, outer = TRUE)

dev.off()

#endregion

#region ############# PLOTS WITH DIFFERENT FILTERS - GOSAT S #############

### STATS ###

# Mean up the GOSAT polarizations
mean_GOSAT_SIF_740nm_unfiltered          <- unfiltered_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_veg                 <- veg_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_time_veg            <- time_veg_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_vza_time_veg        <- vza_time_veg_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_n_vza_time_veg      <- n_vza_time_veg_matched_sounding_means$SIF_740nm_S
mean_GOSAT_SIF_740nm_temp_n_vza_time_veg <- temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_S

# Variances of the errors
var_err_GOSAT_SIF_740nm_unfiltered          <- var(unfiltered_matched_sounding_means$SIF_Uncertainty_740nm_S)
var_err_GOSAT_SIF_740nm_veg                 <- var(veg_matched_sounding_means$SIF_Uncertainty_740nm_S)
var_err_GOSAT_SIF_740nm_time_veg            <- var(time_veg_matched_sounding_means$SIF_Uncertainty_740nm_S)
var_err_GOSAT_SIF_740nm_vza_time_veg        <- var(vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_S)
var_err_GOSAT_SIF_740nm_n_vza_time_veg      <- var(n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_S)
var_err_GOSAT_SIF_740nm_temp_n_vza_time_veg <- var(temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_S)

var_err_OCO2_SIF_740nm_unfiltered          <- var(unfiltered_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_veg                 <- var(veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_time_veg            <- var(time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_vza_time_veg        <- var(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_n_vza_time_veg      <- var(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_temp_n_vza_time_veg <- var(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)

var_err_ratio_SIF_740nm_unfiltered          <- var_err_OCO2_SIF_740nm_unfiltered / var_err_GOSAT_SIF_740nm_unfiltered
var_err_ratio_SIF_740nm_veg                 <- var_err_OCO2_SIF_740nm_veg / var_err_GOSAT_SIF_740nm_veg
var_err_ratio_SIF_740nm_time_veg            <- var_err_OCO2_SIF_740nm_time_veg / var_err_GOSAT_SIF_740nm_time_veg 
var_err_ratio_SIF_740nm_vza_time_veg        <- var_err_OCO2_SIF_740nm_vza_time_veg / var_err_GOSAT_SIF_740nm_vza_time_veg
var_err_ratio_SIF_740nm_n_vza_time_veg      <- var_err_OCO2_SIF_740nm_n_vza_time_veg / var_err_GOSAT_SIF_740nm_n_vza_time_veg
var_err_ratio_SIF_740nm_temp_n_vza_time_veg <- var_err_OCO2_SIF_740nm_temp_n_vza_time_veg / var_err_GOSAT_SIF_740nm_temp_n_vza_time_veg

# Run Deming regressions
reg_SIF_740nm_unfiltered          <- mcreg(mean_GOSAT_SIF_740nm_unfiltered, unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_unfiltered,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_veg                 <- mcreg(mean_GOSAT_SIF_740nm_veg, veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_time_veg            <- mcreg(mean_GOSAT_SIF_740nm_time_veg, time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_vza_time_veg        <- mcreg(mean_GOSAT_SIF_740nm_vza_time_veg, vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_vza_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_n_vza_time_veg      <- mcreg(mean_GOSAT_SIF_740nm_n_vza_time_veg, n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_n_vza_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_temp_n_vza_time_veg <- mcreg(mean_GOSAT_SIF_740nm_temp_n_vza_time_veg, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_temp_n_vza_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")

# Run linear regressions for R2 values
lm_SIF_740nm_unfiltered          <- lm(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered)
lm_SIF_740nm_veg                 <- lm(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg)
lm_SIF_740nm_time_veg            <- lm(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg)
lm_SIF_740nm_vza_time_veg        <- lm(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg)
lm_SIF_740nm_n_vza_time_veg      <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg)
lm_SIF_740nm_temp_n_vza_time_veg <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg)

# Summaries
summary_unfiltered          <- summary(lm_SIF_740nm_unfiltered)
summary_veg                 <- summary(lm_SIF_740nm_veg)
summary_time_veg            <- summary(lm_SIF_740nm_time_veg)
summary_vza_time_veg        <- summary(lm_SIF_740nm_vza_time_veg)
summary_n_vza_time_veg      <- summary(lm_SIF_740nm_n_vza_time_veg)
summary_temp_n_vza_time_veg <- summary(lm_SIF_740nm_temp_n_vza_time_veg)

# Round coefficients
cf_lm_SIF_740nm_unfiltered           <- round(coef(lm_SIF_740nm_unfiltered), 2)
cf_lm_SIF_740nm_veg                  <- round(coef(lm_SIF_740nm_veg), 2)
cf_lm_SIF_740nm_time_veg             <- round(coef(lm_SIF_740nm_time_veg), 2)
cf_lm_SIF_740nm_vza_time_veg         <- round(coef(lm_SIF_740nm_vza_time_veg), 2)
cf_lm_SIF_740nm_n_vza_time_veg       <- round(coef(lm_SIF_740nm_n_vza_time_veg), 2)
cf_lm_SIF_740nm_temp_n_vza_time_veg  <- round(coef(lm_SIF_740nm_temp_n_vza_time_veg), 2)

# Equations
eq_unfiltered          <- paste0("y = ", round(reg_SIF_740nm_unfiltered@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_unfiltered@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_unfiltered@para[2], 2)), "x")
eq_veg                 <- paste0("y = ", round(reg_SIF_740nm_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_veg@para[2], 2)), "x")
eq_time_veg            <- paste0("y = ", round(reg_SIF_740nm_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_time_veg@para[2], 2)), "x")
eq_vza_time_veg        <- paste0("y = ", round(reg_SIF_740nm_vza_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_vza_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_vza_time_veg@para[2], 2)), "x")
eq_n_vza_time_veg      <- paste0("y = ", round(reg_SIF_740nm_n_vza_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_n_vza_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_n_vza_time_veg@para[2], 2)), "x")
eq_temp_n_vza_time_veg <- paste0("y = ", round(reg_SIF_740nm_temp_n_vza_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_temp_n_vza_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_temp_n_vza_time_veg@para[2], 2)), "x")

# R values
r_unfiltered          <- paste0("p < .001, R2 = ", round(summary_unfiltered$adj.r.squared, 2))
r_veg                 <- paste0("p < .001, R2 = ", round(summary_veg$adj.r.squared, 2))
r_time_veg            <- paste0("p < .001, R2 = ", round(summary_time_veg$adj.r.squared, 2))
r_vza_time_veg        <- paste0("p < .001, R2 = ", round(summary_vza_time_veg$adj.r.squared, 2))
r_n_vza_time_veg      <- paste0("p < .001, R2 = ", round(summary_n_vza_time_veg$adj.r.squared, 2))
r_temp_n_vza_time_veg <- paste0("p < .001, R2 = ", round(summary_temp_n_vza_time_veg$adj.r.squared, 2))


##### SETUP MAIN PLOT #####

pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_Deming_GOSAT_S_OCO2_SIF740.pdf", width=8.5, height=6, compress=FALSE)

par(mfrow=c(2, 3), oma = c(0, 1.5, 0, 3.5))

# FILTERED DATA PLOTS

# Unfiltered
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_unfiltered, x.lab = "", y.lab = "", main = title_1, xlim = c(-6, 6), ylim = c(-6, 6),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_unfiltered, r_unfiltered), bty = 'o', bg = "white", box.col = "white")

# Veg
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_veg, x.lab = "", y.lab = "", main = title_2, xlim = c(-6, 6), ylim = c(-6, 6),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_veg, r_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_time_veg, x.lab = "", y.lab = "", main = title_3, xlim = c(-6, 6), ylim = c(-6, 6),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_time_veg, r_time_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time VZA
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_vza_time_veg, x.lab = "", y.lab = "", main = title_4, xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_vza_time_veg, r_vza_time_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time VZA N
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_n_vza_time_veg, x.lab = "", y.lab = "", main = title_5, xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_n_vza_time_veg, r_n_vza_time_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time VZA N Temp
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_temp_n_vza_time_veg, x.lab = "", y.lab = "", main = title_6, xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_temp_n_vza_time_veg, r_temp_n_vza_time_veg), bty = 'o', bg = "white", box.col = "white")

mtext(expression(paste("Mean OCO SIF 740nm")), 2, -0.5, cex = 1, outer = TRUE)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, -1.5, cex = 1, outer = TRUE)

dev.off()


#endregion

#region ############# PLOTS WITH DIFFERENT FILTERS - GOSAT Mean PS #############

### STATS ###

# Mean up the GOSAT polarizations
mean_GOSAT_SIF_740nm_unfiltered          <- (unfiltered_matched_sounding_means$SIF_740nm_P + unfiltered_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_veg                 <- (veg_matched_sounding_means$SIF_740nm_P + veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_time_veg            <- (time_veg_matched_sounding_means$SIF_740nm_P + time_veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_vza_time_veg        <- (vza_time_veg_matched_sounding_means$SIF_740nm_P + vza_time_veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_n_vza_time_veg      <- (n_vza_time_veg_matched_sounding_means$SIF_740nm_P + n_vza_time_veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_740nm_temp_n_vza_time_veg <- (temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_S) / 2


# Variances of the errors
var_err_GOSAT_SIF_740nm_unfiltered          <- var((unfiltered_matched_sounding_means$SIF_Uncertainty_740nm_P + unfiltered_matched_sounding_means$SIF_Uncertainty_740nm_S) / 2)
var_err_GOSAT_SIF_740nm_veg                 <- var((veg_matched_sounding_means$SIF_Uncertainty_740nm_P + veg_matched_sounding_means$SIF_Uncertainty_740nm_S) / 2)
var_err_GOSAT_SIF_740nm_time_veg            <- var((time_veg_matched_sounding_means$SIF_Uncertainty_740nm_P + time_veg_matched_sounding_means$SIF_Uncertainty_740nm_S) / 2)
var_err_GOSAT_SIF_740nm_vza_time_veg        <- var((vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_P + vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_S) / 2)
var_err_GOSAT_SIF_740nm_n_vza_time_veg      <- var((n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_P + n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_S) / 2)
var_err_GOSAT_SIF_740nm_temp_n_vza_time_veg <- var((temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_S) / 2)

var_err_OCO2_SIF_740nm_unfiltered          <- var(unfiltered_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_veg                 <- var(veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_time_veg            <- var(time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_vza_time_veg        <- var(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_n_vza_time_veg      <- var(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_740nm_temp_n_vza_time_veg <- var(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)

var_err_ratio_SIF_740nm_unfiltered          <- var_err_OCO2_SIF_740nm_unfiltered / var_err_GOSAT_SIF_740nm_unfiltered
var_err_ratio_SIF_740nm_veg                 <- var_err_OCO2_SIF_740nm_veg / var_err_GOSAT_SIF_740nm_veg
var_err_ratio_SIF_740nm_time_veg            <- var_err_OCO2_SIF_740nm_time_veg / var_err_GOSAT_SIF_740nm_time_veg 
var_err_ratio_SIF_740nm_vza_time_veg        <- var_err_OCO2_SIF_740nm_vza_time_veg / var_err_GOSAT_SIF_740nm_vza_time_veg
var_err_ratio_SIF_740nm_n_vza_time_veg      <- var_err_OCO2_SIF_740nm_n_vza_time_veg / var_err_GOSAT_SIF_740nm_n_vza_time_veg
var_err_ratio_SIF_740nm_temp_n_vza_time_veg <- var_err_OCO2_SIF_740nm_temp_n_vza_time_veg / var_err_GOSAT_SIF_740nm_temp_n_vza_time_veg

# Run Deming regressions
reg_SIF_740nm_unfiltered          <- mcreg(mean_GOSAT_SIF_740nm_unfiltered, unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_unfiltered,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_veg                 <- mcreg(mean_GOSAT_SIF_740nm_veg, veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_time_veg            <- mcreg(mean_GOSAT_SIF_740nm_time_veg, time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_vza_time_veg        <- mcreg(mean_GOSAT_SIF_740nm_vza_time_veg, vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_vza_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_n_vza_time_veg      <- mcreg(mean_GOSAT_SIF_740nm_n_vza_time_veg, n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_n_vza_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
reg_SIF_740nm_temp_n_vza_time_veg <- mcreg(mean_GOSAT_SIF_740nm_temp_n_vza_time_veg, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm_temp_n_vza_time_veg,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")

# Run linear regressions for R2 values
lm_SIF_740nm_unfiltered          <- lm(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered)
lm_SIF_740nm_veg                 <- lm(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg)
lm_SIF_740nm_time_veg            <- lm(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg)
lm_SIF_740nm_vza_time_veg        <- lm(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg)
lm_SIF_740nm_n_vza_time_veg      <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg)
lm_SIF_740nm_temp_n_vza_time_veg <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg)

# Summaries
summary_unfiltered          <- summary(lm_SIF_740nm_unfiltered)
summary_veg                 <- summary(lm_SIF_740nm_veg)
summary_time_veg            <- summary(lm_SIF_740nm_time_veg)
summary_vza_time_veg        <- summary(lm_SIF_740nm_vza_time_veg)
summary_n_vza_time_veg      <- summary(lm_SIF_740nm_n_vza_time_veg)
summary_temp_n_vza_time_veg <- summary(lm_SIF_740nm_temp_n_vza_time_veg)

# Round coefficients
cf_lm_SIF_740nm_unfiltered           <- round(coef(lm_SIF_740nm_unfiltered), 2)
cf_lm_SIF_740nm_veg                  <- round(coef(lm_SIF_740nm_veg), 2)
cf_lm_SIF_740nm_time_veg             <- round(coef(lm_SIF_740nm_time_veg), 2)
cf_lm_SIF_740nm_vza_time_veg         <- round(coef(lm_SIF_740nm_vza_time_veg), 2)
cf_lm_SIF_740nm_n_vza_time_veg       <- round(coef(lm_SIF_740nm_n_vza_time_veg), 2)
cf_lm_SIF_740nm_temp_n_vza_time_veg  <- round(coef(lm_SIF_740nm_temp_n_vza_time_veg), 2)

# Equations
eq_unfiltered          <- paste0("y = ", round(reg_SIF_740nm_unfiltered@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_unfiltered@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_unfiltered@para[2], 2)), "x")
eq_veg                 <- paste0("y = ", round(reg_SIF_740nm_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_veg@para[2], 2)), "x")
eq_time_veg            <- paste0("y = ", round(reg_SIF_740nm_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_time_veg@para[2], 2)), "x")
eq_vza_time_veg        <- paste0("y = ", round(reg_SIF_740nm_vza_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_vza_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_vza_time_veg@para[2], 2)), "x")
eq_n_vza_time_veg      <- paste0("y = ", round(reg_SIF_740nm_n_vza_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_n_vza_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_n_vza_time_veg@para[2], 2)), "x")
eq_temp_n_vza_time_veg <- paste0("y = ", round(reg_SIF_740nm_temp_n_vza_time_veg@para[1], 2),
                                 ifelse(sign(round(reg_SIF_740nm_temp_n_vza_time_veg@para[2], 2)) == 1, " + ", " - "), abs(round(reg_SIF_740nm_temp_n_vza_time_veg@para[2], 2)), "x")

# R values
r_unfiltered          <- paste0("p < .001, R2 = ", round(summary_unfiltered$adj.r.squared, 2))
r_veg                 <- paste0("p < .001, R2 = ", round(summary_veg$adj.r.squared, 2))
r_time_veg            <- paste0("p < .001, R2 = ", round(summary_time_veg$adj.r.squared, 2))
r_vza_time_veg        <- paste0("p < .001, R2 = ", round(summary_vza_time_veg$adj.r.squared, 2))
r_n_vza_time_veg      <- paste0("p < .001, R2 = ", round(summary_n_vza_time_veg$adj.r.squared, 2))
r_temp_n_vza_time_veg <- paste0("p < .001, R2 = ", round(summary_temp_n_vza_time_veg$adj.r.squared, 2))


##### SETUP MAIN PLOT #####

pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_Deming_GOSAT_PS_OCO2_SIF740.pdf", width=8.5, height=6, compress=FALSE)

par(mfrow=c(2, 3), oma = c(0, 1.5, 0, 3.5))

# FILTERED DATA PLOTS

# Unfiltered
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_unfiltered, x.lab = "", y.lab = "", main = title_1, xlim = c(-6, 6), ylim = c(-6, 6),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_unfiltered, r_unfiltered), bty = 'o', bg = "white", box.col = "white")

# Veg
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_veg, x.lab = "", y.lab = "", main = title_2, xlim = c(-6, 6), ylim = c(-6, 6),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_veg, r_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_time_veg, x.lab = "", y.lab = "", main = title_3, xlim = c(-6, 6), ylim = c(-6, 6),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_time_veg, r_time_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time VZA
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_vza_time_veg, x.lab = "", y.lab = "", main = title_4, xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_vza_time_veg, r_vza_time_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time VZA N
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_n_vza_time_veg, x.lab = "", y.lab = "", main = title_5, xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_n_vza_time_veg, r_n_vza_time_veg), bty = 'o', bg = "white", box.col = "white")

# Veg Time VZA N Temp
op <- par(mar = c(3.5, 2, 2, 0))
MCResult.plot(reg_SIF_740nm_temp_n_vza_time_veg, x.lab = "", y.lab = "", main = title_6, xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_temp_n_vza_time_veg, r_temp_n_vza_time_veg), bty = 'o', bg = "white", box.col = "white")

mtext(expression(paste("Mean OCO SIF 740nm")), 2, -0.5, cex = 1, outer = TRUE)
mtext(expression(paste("Mean GOSAT PS SIF 740nm")), 1, -1.5, cex = 1, outer = TRUE)

dev.off()

#endregion


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


##### HISTOGRAM OF TIME ######

hist(veg_matched_sounding_means$Delta_Time, breaks = "quarters", freq = TRUE, format = "%Y-%m", right = FALSE)
