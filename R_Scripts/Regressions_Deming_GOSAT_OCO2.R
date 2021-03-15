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

# Run linear regressions
reg_SIF_740nm_unfiltered          <- lm(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered)
reg_SIF_740nm_veg                 <- lm(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg)
reg_SIF_740nm_time_veg            <- lm(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg)
reg_SIF_740nm_vza_time_veg        <- lm(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg)
reg_SIF_740nm_n_vza_time_veg      <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg)
reg_SIF_740nm_temp_n_vza_time_veg <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg)


# Summaries
summary_unfiltered          <- summary(reg_SIF_740nm_unfiltered)
summary_veg                 <- summary(reg_SIF_740nm_veg)
summary_time_veg            <- summary(reg_SIF_740nm_time_veg)
summary_vza_time_veg        <- summary(reg_SIF_740nm_vza_time_veg)
summary_n_vza_time_veg      <- summary(reg_SIF_740nm_n_vza_time_veg)
summary_temp_n_vza_time_veg <- summary(reg_SIF_740nm_temp_n_vza_time_veg)

# Round coefficients
cf_reg_SIF_740nm_unfiltered           <- round(coef(reg_SIF_740nm_unfiltered), 2)
cf_reg_SIF_740nm_veg                  <- round(coef(reg_SIF_740nm_veg), 2)
cf_reg_SIF_740nm_time_veg             <- round(coef(reg_SIF_740nm_time_veg), 2)
cf_reg_SIF_740nm_vza_time_veg         <- round(coef(reg_SIF_740nm_vza_time_veg), 2)
cf_reg_SIF_740nm_n_vza_time_veg       <- round(coef(reg_SIF_740nm_n_vza_time_veg), 2)
cf_reg_SIF_740nm_temp_n_vza_time_veg  <- round(coef(reg_SIF_740nm_temp_n_vza_time_veg), 2)

# Equations
eq_unfiltered     <- paste0("y = ", cf_reg_SIF_740nm_unfiltered[1],
                            ifelse(sign(cf_reg_SIF_740nm_unfiltered[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_unfiltered[2]), "x")
eq_veg            <- paste0("y = ", cf_reg_SIF_740nm_veg[1],
                            ifelse(sign(cf_reg_SIF_740nm_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_veg[2]), "x")
eq_time_veg       <- paste0("y = ", cf_reg_SIF_740nm_time_veg[1],
                            ifelse(sign(cf_reg_SIF_740nm_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_time_veg[2]), "x")
eq_vza_time_veg   <- paste0("y = ", cf_reg_SIF_740nm_vza_time_veg[1],
                            ifelse(sign(cf_reg_SIF_740nm_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_vza_time_veg[2]), "x")
eq_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_n_vza_time_veg[1],
                            ifelse(sign(cf_reg_SIF_740nm_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_n_vza_time_veg[2]), "x")
eq_temp_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_temp_n_vza_time_veg[1],
                                 ifelse(sign(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]), "x")

# RMSE
RMSE_unfiltered          <- paste0("RMSE = ", round(sqrt(mean(summary_unfiltered$residuals ^ 2)), 2))
RMSE_veg                 <- paste0("RMSE = ", round(sqrt(mean(summary_veg$residuals ^ 2)), 2))
RMSE_time_veg            <- paste0("RMSE = ", round(sqrt(mean(summary_time_veg$residuals ^ 2)), 2))
RMSE_vza_time_veg        <- paste0("RMSE = ", round(sqrt(mean(summary_vza_time_veg$residuals ^ 2)), 2))
RMSE_n_vza_time_veg      <- paste0("RMSE = ", round(sqrt(mean(summary_n_vza_time_veg$residuals ^ 2)), 2))
RMSE_temp_n_vza_time_veg <- paste0("RMSE = ", round(sqrt(mean(summary_temp_n_vza_time_veg$residuals ^ 2)), 2))

# R values
r_unfiltered          <- paste0("p < .001, R2 = ", round(summary_unfiltered$adj.r.squared, 2))
r_veg                 <- paste0("p < .001, R2 = ", round(summary_veg$adj.r.squared, 2))
r_time_veg            <- paste0("p < .001, R2 = ", round(summary_time_veg$adj.r.squared, 2))
r_vza_time_veg        <- paste0("p < .001, R2 = ", round(summary_vza_time_veg$adj.r.squared, 2))
r_n_vza_time_veg      <- paste0("p < .001, R2 = ", round(summary_n_vza_time_veg$adj.r.squared, 2))
r_temp_n_vza_time_veg <- paste0("p < .001, R2 = ", round(summary_temp_n_vza_time_veg$adj.r.squared, 2))

##### SETUP MAIN PLOT #####

# pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_GOSAT_S_OCO2_SIF740.pdf", width=8, height=8, compress=FALSE)

par(mfrow=c(2, 3), oma = c(0, 0.5, 0, 1.5))

# FILTERED DATA PLOTS

# Unfiltered
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_unfiltered, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_unfiltered, r_unfiltered, RMSE_unfiltered), bty = 'n')
mtext(expression(paste("Unfiltered")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_veg, r_veg, RMSE_veg), bty = 'n')
mtext(expression(paste("Vegetation Only")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_time_veg, r_time_veg, RMSE_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_vza_time_veg, r_vza_time_veg, RMSE_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_n_vza_time_veg, r_n_vza_time_veg, RMSE_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N Temp
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_temp_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_temp_n_vza_time_veg, r_temp_n_vza_time_veg, RMSE_temp_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10 / T >= 5C")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("GOSAT S SIF 740nm")), 1, 1.5, cex = 1)
box()

# dev.off()

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

# Run linear regressions
reg_SIF_740nm_unfiltered          <- lm(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered)
reg_SIF_740nm_veg                 <- lm(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg)
reg_SIF_740nm_time_veg            <- lm(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg)
reg_SIF_740nm_vza_time_veg        <- lm(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg)
reg_SIF_740nm_n_vza_time_veg      <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg)
reg_SIF_740nm_temp_n_vza_time_veg <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg)


# Summaries
summary_unfiltered          <- summary(reg_SIF_740nm_unfiltered)
summary_veg                 <- summary(reg_SIF_740nm_veg)
summary_time_veg            <- summary(reg_SIF_740nm_time_veg)
summary_vza_time_veg        <- summary(reg_SIF_740nm_vza_time_veg)
summary_n_vza_time_veg      <- summary(reg_SIF_740nm_n_vza_time_veg)
summary_temp_n_vza_time_veg <- summary(reg_SIF_740nm_temp_n_vza_time_veg)

# Round coefficients
cf_reg_SIF_740nm_unfiltered           <- round(coef(reg_SIF_740nm_unfiltered), 2)
cf_reg_SIF_740nm_veg                  <- round(coef(reg_SIF_740nm_veg), 2)
cf_reg_SIF_740nm_time_veg             <- round(coef(reg_SIF_740nm_time_veg), 2)
cf_reg_SIF_740nm_vza_time_veg         <- round(coef(reg_SIF_740nm_vza_time_veg), 2)
cf_reg_SIF_740nm_n_vza_time_veg       <- round(coef(reg_SIF_740nm_n_vza_time_veg), 2)
cf_reg_SIF_740nm_temp_n_vza_time_veg  <- round(coef(reg_SIF_740nm_temp_n_vza_time_veg), 2)

# Equations
eq_unfiltered     <- paste0("y = ", cf_reg_SIF_740nm_unfiltered[1],
                            ifelse(sign(cf_reg_SIF_740nm_unfiltered[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_unfiltered[2]), "x")
eq_veg            <- paste0("y = ", cf_reg_SIF_740nm_veg[1],
                            ifelse(sign(cf_reg_SIF_740nm_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_veg[2]), "x")
eq_time_veg       <- paste0("y = ", cf_reg_SIF_740nm_time_veg[1],
                            ifelse(sign(cf_reg_SIF_740nm_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_time_veg[2]), "x")
eq_vza_time_veg   <- paste0("y = ", cf_reg_SIF_740nm_vza_time_veg[1],
                            ifelse(sign(cf_reg_SIF_740nm_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_vza_time_veg[2]), "x")
eq_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_n_vza_time_veg[1],
                            ifelse(sign(cf_reg_SIF_740nm_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_n_vza_time_veg[2]), "x")
eq_temp_n_vza_time_veg <- paste0("y = ", cf_reg_SIF_740nm_temp_n_vza_time_veg[1],
                                 ifelse(sign(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm_temp_n_vza_time_veg[2]), "x")

# RMSE
RMSE_unfiltered          <- paste0("RMSE = ", round(sqrt(mean(summary_unfiltered$residuals ^ 2)), 2))
RMSE_veg                 <- paste0("RMSE = ", round(sqrt(mean(summary_veg$residuals ^ 2)), 2))
RMSE_time_veg            <- paste0("RMSE = ", round(sqrt(mean(summary_time_veg$residuals ^ 2)), 2))
RMSE_vza_time_veg        <- paste0("RMSE = ", round(sqrt(mean(summary_vza_time_veg$residuals ^ 2)), 2))
RMSE_n_vza_time_veg      <- paste0("RMSE = ", round(sqrt(mean(summary_n_vza_time_veg$residuals ^ 2)), 2))
RMSE_temp_n_vza_time_veg <- paste0("RMSE = ", round(sqrt(mean(summary_temp_n_vza_time_veg$residuals ^ 2)), 2))

# R values
r_unfiltered          <- paste0("p < .001, R2 = ", round(summary_unfiltered$adj.r.squared, 2))
r_veg                 <- paste0("p < .001, R2 = ", round(summary_veg$adj.r.squared, 2))
r_time_veg            <- paste0("p < .001, R2 = ", round(summary_time_veg$adj.r.squared, 2))
r_vza_time_veg        <- paste0("p < .001, R2 = ", round(summary_vza_time_veg$adj.r.squared, 2))
r_n_vza_time_veg      <- paste0("p < .001, R2 = ", round(summary_n_vza_time_veg$adj.r.squared, 2))
r_temp_n_vza_time_veg <- paste0("p < .001, R2 = ", round(summary_temp_n_vza_time_veg$adj.r.squared, 2))

##### SETUP MAIN PLOT #####

# pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_GOSAT_OCO2_SIF740.pdf", width=8, height=8, compress=FALSE)

par(mfrow=c(2, 3), oma = c(0, 0.5, 0, 1.5))

# FILTERED DATA PLOTS

# Unfiltered
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(unfiltered_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_unfiltered, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_unfiltered, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_unfiltered, r_unfiltered, RMSE_unfiltered), bty = 'n')
mtext(expression(paste("Unfiltered")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_veg, r_veg, RMSE_veg), bty = 'n')
mtext(expression(paste("Vegetation Only")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_time_veg, r_time_veg, RMSE_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_vza_time_veg, r_vza_time_veg, RMSE_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_n_vza_time_veg, r_n_vza_time_veg, RMSE_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# Veg Time VZA N Temp
op <- par(mar = c(3.5, 3.5, 2, 0))
reg_plot <- plot(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm_temp_n_vza_time_veg, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm_temp_n_vza_time_veg, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_temp_n_vza_time_veg, r_temp_n_vza_time_veg, RMSE_temp_n_vza_time_veg), bty = 'n')
mtext(expression(paste("Veg / < 1 hour / VZA < 5 / N >= 10 / T >= 5C")), 3, 0.5, cex = 1)
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# dev.off()

#endregion

#region ############# SIF, DAILY, RELATIVE ################

### STATS ###

# Mean up the GOSAT polarizations
mean_GOSAT_SIF_740nm          <- (n_vza_time_veg_matched_sounding_means$SIF_740nm_P + n_vza_time_veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_771nm          <- (n_vza_time_veg_matched_sounding_means$SIF_757nm_P + n_vza_time_veg_matched_sounding_means$SIF_757nm_S) / 2
mean_GOSAT_SIF_771nm          <- (n_vza_time_veg_matched_sounding_means$SIF_771nm_P + n_vza_time_veg_matched_sounding_means$SIF_771nm_S) / 2
mean_GOSAT_SIF_Relative_757nm <- (n_vza_time_veg_matched_sounding_means$SIF_Relative_757nm_P + n_vza_time_veg_matched_sounding_means$SIF_Relative_757nm_S) / 2
mean_GOSAT_SIF_Relative_771nm <- (n_vza_time_veg_matched_sounding_means$SIF_Relative_771nm_P + n_vza_time_veg_matched_sounding_means$SIF_Relative_771nm_S) / 2
mean_GOSAT_SIF_Daily_740nm    <- (n_vza_time_veg_matched_sounding_means$SIF_Daily_740nm_P + n_vza_time_veg_matched_sounding_means$SIF_Daily_740nm_S) / 2
mean_GOSAT_SIF_Daily_757nm    <- (n_vza_time_veg_matched_sounding_means$SIF_Daily_757nm_P + n_vza_time_veg_matched_sounding_means$SIF_Daily_757nm_S) / 2
mean_GOSAT_SIF_Daily_771nm    <- (n_vza_time_veg_matched_sounding_means$SIF_Daily_771nm_P + n_vza_time_veg_matched_sounding_means$SIF_Daily_771nm_S) / 2

# Run linear regressions
reg_SIF_740nm          <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm)
reg_SIF_757nm          <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_757nm ~ mean_GOSAT_SIF_757nm)
reg_SIF_771nm          <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_771nm ~ mean_GOSAT_SIF_771nm)
reg_SIF_Relative_757nm <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_757nm ~ mean_GOSAT_SIF_Relative_757nm)
reg_SIF_Relative_771nm <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_771nm ~ mean_GOSAT_SIF_Relative_771nm)
reg_SIF_Daily_740nm    <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_740nm ~ mean_GOSAT_SIF_Daily_740nm)
reg_SIF_Daily_757nm    <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_757nm ~ mean_GOSAT_SIF_Daily_757nm)
reg_SIF_Daily_771nm    <- lm(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_771nm ~ mean_GOSAT_SIF_Daily_771nm)

# Summaries
summary_SIF_740nm          <- summary(reg_SIF_740nm)
summary_SIF_757nm          <- summary(reg_SIF_757nm)
summary_SIF_771nm          <- summary(reg_SIF_771nm)
summary_SIF_Relative_757nm <- summary(reg_SIF_Relative_757nm)
summary_SIF_Relative_771nm <- summary(reg_SIF_Relative_771nm)
summary_SIF_Daily_740nm    <- summary(reg_SIF_Daily_740nm)
summary_SIF_Daily_757nm    <- summary(reg_SIF_Daily_757nm)
summary_SIF_Daily_771nm    <- summary(reg_SIF_Daily_771nm)

# Round coefficients
cf_reg_SIF_740nm <- round(coef(reg_SIF_740nm), 2)
cf_reg_SIF_757nm <- round(coef(reg_SIF_757nm), 2)
cf_reg_SIF_771nm <- round(coef(reg_SIF_771nm), 2)
cf_reg_SIF_Relative_757nm <- round(coef(reg_SIF_Relative_757nm), 2)
cf_reg_SIF_Relative_771nm <- round(coef(reg_SIF_Relative_771nm), 2)
cf_reg_SIF_Daily_740nm <- round(coef(reg_SIF_Daily_740nm), 2)
cf_reg_SIF_Daily_757nm <- round(coef(reg_SIF_Daily_757nm), 2)
cf_reg_SIF_Daily_771nm <- round(coef(reg_SIF_Daily_771nm), 2)

# Equations
eq_SIF_740nm          <- paste0("y = ", cf_reg_SIF_740nm[1],
                                ifelse(sign(cf_reg_SIF_740nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_740nm[2]), "x")
eq_SIF_757nm          <- paste0("y = ", cf_reg_SIF_757nm[1],
                                ifelse(sign(cf_reg_SIF_757nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_757nm[2]), "x")
eq_SIF_771nm          <- paste0("y = ", cf_reg_SIF_771nm[1],
                                ifelse(sign(cf_reg_SIF_771nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_771nm[2]), "x")
eq_SIF_Relative_757nm <- paste0("y = ", cf_reg_SIF_Relative_757nm[1],
                                ifelse(sign(cf_reg_SIF_757nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Relative_757nm[2]), "x")
eq_SIF_Relative_771nm <- paste0("y = ", cf_reg_SIF_Relative_771nm[1],
                                ifelse(sign(cf_reg_SIF_771nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Relative_771nm[2]), "x")
eq_SIF_Daily_740nm    <- paste0("y = ", cf_reg_SIF_Daily_740nm[1],
                                ifelse(sign(cf_reg_SIF_Daily_740nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Daily_740nm[2]), "x")
eq_SIF_Daily_757nm    <- paste0("y = ", cf_reg_SIF_Daily_757nm[1],
                                ifelse(sign(cf_reg_SIF_Daily_757nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Daily_757nm[2]), "x")
eq_SIF_Daily_771nm    <- paste0("y = ", cf_reg_SIF_Daily_771nm[1],
                                ifelse(sign(cf_reg_SIF_Daily_771nm[2]) == 1, " + ", " - "), abs(cf_reg_SIF_Daily_771nm[2]), "x")

# RMSE
RMSE_SIF_740nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_740nm$residuals ^ 2)), 2))
RMSE_SIF_757nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_757nm$residuals ^ 2)), 2))
RMSE_SIF_771nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_771nm$residuals ^ 2)), 2))
RMSE_SIF_Relative_757nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Relative_757nm$residuals ^ 2)), 2))
RMSE_SIF_Relative_771nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Relative_771nm$residuals ^ 2)), 2))
RMSE_SIF_Daily_740nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Daily_740nm$residuals ^ 2)), 2))
RMSE_SIF_Daily_757nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Daily_757nm$residuals ^ 2)), 2))
RMSE_SIF_Daily_771nm <- paste0("RMSE = ", round(sqrt(mean(summary_SIF_Daily_771nm$residuals ^ 2)), 2))

# R values
r_SIF_740nm <- paste0("p < .001, R2 = ", round(summary_SIF_740nm$adj.r.squared, 2))
r_SIF_757nm <- paste0("p < .001, R2 = ", round(summary_SIF_757nm$adj.r.squared, 2))
r_SIF_771nm <- paste0("p < .001, R2 = ", round(summary_SIF_771nm$adj.r.squared, 2))
r_SIF_Relative_757nm <- paste0("p < .001, R2 = ", round(summary_SIF_Relative_757nm$adj.r.squared, 2))
r_SIF_Relative_771nm <- paste0("p < .001, R2 = ", round(summary_SIF_Relative_771nm$adj.r.squared, 2))
r_SIF_Daily_740nm <- paste0("p < .001, R2 = ", round(summary_SIF_Daily_740nm$adj.r.squared, 2))
r_SIF_Daily_757nm <- paste0("p < .001, R2 = ", round(summary_SIF_Daily_757nm$adj.r.squared, 2))
r_SIF_Daily_771nm <- paste0("p < .001, R2 = ", round(summary_SIF_Daily_771nm$adj.r.squared, 2))

##### SETUP MAIN PLOT #####

# pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_GOSAT_OCO2.pdf", width=8, height=8, compress=FALSE)

par(mfrow=c(3, 3), oma = c(0, 0.5, 0.5, 0.5))

# FILTERED DATA PLOTS

# SIF 757
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_757nm ~ mean_GOSAT_SIF_757nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_757nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_757nm, r_SIF_757nm, RMSE_SIF_757nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF 757nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 757nm")), 1, 1.5, cex = 1)
box()

# SIF 771
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_771nm ~ mean_GOSAT_SIF_771nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_771nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_771nm, r_SIF_771nm, RMSE_SIF_771nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF 771nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 771nm")), 1, 1.5, cex = 1)
box()

# SIF 740
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_740nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_740nm, r_SIF_740nm, RMSE_SIF_740nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF 740nm")), 1, 1.5, cex = 1)
box()

# SIF DAILY 757
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_757nm ~ mean_GOSAT_SIF_Daily_757nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Daily_757nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Daily_757nm, r_SIF_Daily_757nm, RMSE_SIF_Daily_757nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Daily 757nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Daily 757nm")), 1, 1.5, cex = 1)
box()

# SIF DAILY 771
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_771nm ~ mean_GOSAT_SIF_Daily_771nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Daily_771nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Daily_771nm, r_SIF_Daily_771nm, RMSE_SIF_Daily_771nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Daily 771nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Daily 771nm")), 1, 1.5, cex = 1)
box()

# SIF DAILY 740
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_740nm ~ mean_GOSAT_SIF_Daily_740nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Daily_740nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Daily_740nm, r_SIF_Daily_740nm, RMSE_SIF_Daily_740nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Daily 740nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Daily 740nm")), 1, 1.5, cex = 1)
box()

# SIF Relative 757
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_757nm ~ mean_GOSAT_SIF_Relative_757nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Relative_757nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Relative_757nm, r_SIF_Relative_757nm, RMSE_SIF_Relative_757nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Relative 757nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Relative 757nm")), 1, 1.5, cex = 1)
box()

# SIF Relative 771
op <- par(mar = c(3.5, 3.5, 0, 0))
reg_plot <- plot(n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_771nm ~ mean_GOSAT_SIF_Relative_771nm, pch = 1, axes = FALSE, ann = FALSE, xlim = c(-3, 3), ylim = c(-3, 3), cex = 1)
abline(reg_SIF_Relative_771nm, lwd = 1, lty = 1)
axis(side = 1, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
axis(side = 2, tck = 0.03, cex.axis = 1, mgp = c(3, 0.3, 0), at = c(seq(from = -3, to = 3, by = 1)))
legend('topleft', legend = c(eq_SIF_Relative_771nm, r_SIF_Relative_771nm, RMSE_SIF_Relative_771nm), cex = 0.8, bty = 'n')
mtext(expression(paste("Mean OCO SIF Relative 771nm")), 2, 1.5, cex = 1)
mtext(expression(paste("Mean GOSAT SIF Relative 771nm")), 1, 1.5, cex = 1)
box()

# dev.off()
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
