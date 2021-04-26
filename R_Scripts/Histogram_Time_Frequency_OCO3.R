
unfiltered_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/OCO3/unfiltered_matched_sounding_means_oco3.csv")

temp_n_vza_time_veg_matched_sounding_means <- read.csv("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/CSV/OCO3/temp_n_vza_time_veg_matched_sounding_means_oco3.csv")

dates_all    <- as.Date(unfiltered_matched_sounding_means$Delta_Time)
dates_filter <- as.Date(temp_n_vza_time_veg_matched_sounding_means$Delta_Time)

d1 <- as.Date(paste0("201501","01"), "%Y%m%d")
d2 <- as.Date(paste0("202007","01"), "%Y%m%d")

months_all    <- seq(min(dates_all), max(dates_all), by = "month")
months_filter <- seq(min(dates_filter), max(dates_filter), by = "month")


pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Histogram_Time_Frequency_All_Filters.pdf", width=8, height=6, compress=FALSE)

hist(dates_filter, breaks = "month", freq = TRUE, format = "%Y-%m", right = FALSE,
     main = "Monthly Frequency of Overlapping Soundings for GOSAT and OCO-2\n(All Filters Applied)",
     xlab = NA, col = rgb(44,162,95, max = 255), axes = FALSE)
axis(side = 1, at = months_filter[c(1, 13, 25, 37, 49, 61)], labels = format(months_filter[c(1, 13, 25, 37, 49, 61)], "%Y-%m"), mgp = c(1.5, 0.3, 0), tck = 0.01)
axis(side = 2, mgp = c(1.5, 0.3, 0), tck = 0.01)
box()

dev.off()

pdf("C:/Russell/Projects/Geometry/R_Scripts/Figures/Match_GOSAT_OCO/Histogram_Time_Frequency_No_Filters.pdf", width=8, height=6, compress=FALSE)

hist(dates_all, breaks = "month", freq = TRUE, format = "%Y-%m",
     main = "Monthly Frequency of Overlapping Soundings for GOSAT and OCO-2\n(QC = Best; Cloud Free)",
     xlab = NA, col = rgb(67,162,202, max = 255), axes = FALSE)
axis(side = 1, at = months_all[c(1, 13, 25, 37, 49, 61)], labels = format(months_all[c(1, 13, 25, 37, 49, 61)], "%Y-%m"), mgp = c(1.5, 0.3, 0), tck = 0.01)
axis(side = 2, mgp = c(1.5, 0.3, 0), tck = 0.01)
box()

dev.off()
