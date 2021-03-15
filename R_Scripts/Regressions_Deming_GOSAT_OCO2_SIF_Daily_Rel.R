library(mcr)

temp_n_vza_time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/temp_n_vza_time_veg_matched_sounding_means.csv")



#region ############# SIF, DAILY, RELATIVE ################

### STATS ###

# Mean up the GOSAT polarizations
mean_GOSAT_SIF_740nm          <- (temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_740nm_S) / 2
mean_GOSAT_SIF_757nm          <- (temp_n_vza_time_veg_matched_sounding_means$SIF_757nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_757nm_S) / 2
mean_GOSAT_SIF_771nm          <- (temp_n_vza_time_veg_matched_sounding_means$SIF_771nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_771nm_S) / 2
mean_GOSAT_SIF_Relative_757nm <- (temp_n_vza_time_veg_matched_sounding_means$SIF_Relative_757nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_Relative_757nm_S) / 2
mean_GOSAT_SIF_Relative_771nm <- (temp_n_vza_time_veg_matched_sounding_means$SIF_Relative_771nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_Relative_771nm_S) / 2
mean_GOSAT_SIF_Daily_740nm    <- (temp_n_vza_time_veg_matched_sounding_means$SIF_Daily_740nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_Daily_740nm_S) / 2
mean_GOSAT_SIF_Daily_757nm    <- (temp_n_vza_time_veg_matched_sounding_means$SIF_Daily_757nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_Daily_757nm_S) / 2
mean_GOSAT_SIF_Daily_771nm    <- (temp_n_vza_time_veg_matched_sounding_means$SIF_Daily_771nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_Daily_771nm_S) / 2

# Variances of the errors
var_err_GOSAT_SIF_740nm <- var((temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_740nm_S) / 2)
var_err_GOSAT_SIF_757nm <- var((temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_757nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_757nm_S) / 2)
var_err_GOSAT_SIF_771nm <- var((temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_771nm_P + temp_n_vza_time_veg_matched_sounding_means$SIF_Uncertainty_771nm_S) / 2)

var_err_OCO2_SIF_740nm <- var(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_740nm)
var_err_OCO2_SIF_757nm <- var(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_757nm)
var_err_OCO2_SIF_771nm <- var(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Uncertainty_771nm)

var_err_ratio_SIF_740nm <- var_err_OCO2_SIF_740nm / var_err_GOSAT_SIF_740nm
var_err_ratio_SIF_757nm <- var_err_OCO2_SIF_757nm / var_err_GOSAT_SIF_757nm
var_err_ratio_SIF_771nm <- var_err_OCO2_SIF_740nm / var_err_GOSAT_SIF_771nm

# Run Deming regressions
dem_SIF_740nm          <- mcreg(mean_GOSAT_SIF_740nm, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm,
                                           mref.name = "GOSAT", mtest.name = "OCO-2")
dem_SIF_757nm          <- mcreg(mean_GOSAT_SIF_757nm, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_757nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_757nm,
                                mref.name = "GOSAT", mtest.name = "OCO-2")
dem_SIF_771nm          <- mcreg(mean_GOSAT_SIF_771nm, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_771nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_771nm,
                                mref.name = "GOSAT", mtest.name = "OCO-2")
dem_SIF_Relative_757nm          <- mcreg(mean_GOSAT_SIF_Relative_757nm, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_757nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_757nm,
                                mref.name = "GOSAT", mtest.name = "OCO-2")
dem_SIF_Relative_771nm          <- mcreg(mean_GOSAT_SIF_Relative_771nm, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_771nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_771nm,
                                         mref.name = "GOSAT", mtest.name = "OCO-2")
dem_SIF_Daily_740nm          <- mcreg(mean_GOSAT_SIF_Daily_740nm, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_740nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_740nm,
                                mref.name = "GOSAT", mtest.name = "OCO-2")
dem_SIF_Daily_757nm          <- mcreg(mean_GOSAT_SIF_Daily_757nm, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_757nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_757nm,
                                mref.name = "GOSAT", mtest.name = "OCO-2")
dem_SIF_Daily_771nm          <- mcreg(mean_GOSAT_SIF_Daily_771nm, temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_771nm, method.reg = "Deming", error.ratio = var_err_ratio_SIF_771nm,
                                mref.name = "GOSAT", mtest.name = "OCO-2")

# Run linear regressions
reg_SIF_740nm          <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_740nm ~ mean_GOSAT_SIF_740nm)
reg_SIF_757nm          <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_757nm ~ mean_GOSAT_SIF_757nm)
reg_SIF_771nm          <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_771nm ~ mean_GOSAT_SIF_771nm)
reg_SIF_Relative_757nm <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_757nm ~ mean_GOSAT_SIF_Relative_757nm)
reg_SIF_Relative_771nm <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Relative_771nm ~ mean_GOSAT_SIF_Relative_771nm)
reg_SIF_Daily_740nm    <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_740nm ~ mean_GOSAT_SIF_Daily_740nm)
reg_SIF_Daily_757nm    <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_757nm ~ mean_GOSAT_SIF_Daily_757nm)
reg_SIF_Daily_771nm    <- lm(temp_n_vza_time_veg_matched_sounding_means$Mean_OCO_SIF_Daily_771nm ~ mean_GOSAT_SIF_Daily_771nm)

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

pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Regressions_Deming_GOSAT_OCO2_SIF_Daily_Rel.pdf", width=8, height=8, compress=FALSE)

par(mfrow=c(3, 3), oma = c(0, 0.5, 0.5, 0.5))

# FILTERED DATA PLOTS

# SIF #####
op <- par(mar = c(3.5, 3.5, 0, 0))
MCResult.plot(dem_SIF_740nm, x.lab = "Mean GOSAT SIF 740nm", y.lab = "Mean OCO-2 SIF 740nm", main = "", xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_SIF_740nm, r_SIF_740nm), bty = 'o', bg = "white", box.col = "white")

op <- par(mar = c(3.5, 3.5, 0, 0))
MCResult.plot(dem_SIF_757nm, x.lab = "Mean GOSAT SIF 757nm", y.lab = "Mean OCO-2 SIF 757nm", main = "", xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_SIF_757nm, r_SIF_757nm), bty = 'o', bg = "white", box.col = "white")

op <- par(mar = c(3.5, 3.5, 0, 0))
MCResult.plot(dem_SIF_771nm, x.lab = "Mean GOSAT SIF 771nm", y.lab = "Mean OCO-2 SIF 771nm", main = "", xlim = c(-3, 3), ylim = c(-3, 3),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_SIF_771nm, r_SIF_771nm), bty = 'o', bg = "white", box.col = "white")


# SIF Daily #####
op <- par(mar = c(3.5, 3.5, 0, 0))
MCResult.plot(dem_SIF_Daily_740nm, x.lab = "Mean GOSAT SIF Daily 740nm", y.lab = "Mean OCO-2 SIF Daily 740nm", main = "", xlim = c(-1, 1), ylim = c(-1, 1),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_SIF_Daily_740nm, r_SIF_Daily_740nm), bty = 'o', bg = "white", box.col = "white")

op <- par(mar = c(3.5, 3.5, 0, 0))
MCResult.plot(dem_SIF_Daily_757nm, x.lab = "Mean GOSAT SIF Daily 757nm", y.lab = "Mean OCO-2 SIF Daily 757nm", main = "", xlim = c(-1, 1), ylim = c(-1, 1),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_SIF_Daily_757nm, r_SIF_Daily_757nm), bty = 'o', bg = "white", box.col = "white")

op <- par(mar = c(3.5, 3.5, 0, 0))
MCResult.plot(dem_SIF_Daily_771nm, x.lab = "Mean GOSAT SIF Daily 771nm", y.lab = "Mean OCO-2 SIF Daily 771nm", main = "", xlim = c(-1, 1), ylim = c(-1, 1),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_SIF_Daily_771nm, r_SIF_Daily_771nm), bty = 'o', bg = "white", box.col = "white")


# SIF Relative #####
op <- par(mar = c(3.5, 3.5, 0, 0))
MCResult.plot(dem_SIF_Relative_757nm, x.lab = "Mean GOSAT SIF Relative 757nm", y.lab = "Mean OCO-2 SIF Relative 757nm", main = "", xlim = c(-0.025, 0.025), ylim = c(-0.025, 0.025),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_SIF_Relative_757nm, r_SIF_Relative_757nm), bty = 'o', bg = "white", box.col = "white")

op <- par(mar = c(3.5, 3.5, 0, 0))
MCResult.plot(dem_SIF_Relative_771nm, x.lab = "Mean GOSAT SIF Relative 771nm", y.lab = "Mean OCO-2 SIF Relative 771nm", main = "", xlim = c(-0.025, 0.025), ylim = c(-0.025, 0.025),
              sub = "", tck = 0.03, mgp = c(1.5, 0.3, 0), identity.lwd = 2, add.legend = FALSE, add.cor = FALSE)
legend('topleft', inset = c(0.01, 0.01), legend = c(eq_SIF_Relative_771nm, r_SIF_Relative_771nm), bty = 'o', bg = "white", box.col = "white")

dev.off()
#endregion