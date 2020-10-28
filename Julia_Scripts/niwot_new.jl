using CSV
include("Build_OCO3_Data.jl")
include("derive_spectrum.jl")
include("Plot_OCO3_Data.jl")

df_6283 = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_niwot_2020-06-12_6283.csv"
df_6283 = DataFrame!(CSV.File(df_6283));
df_6287 = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_niwot_2020-06-12_6287.csv"
df_6287 = DataFrame!(CSV.File(df_6287));
df_6348 = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_niwot_2020-06-16_6348.csv"
df_6348 = DataFrame!(CSV.File(df_6348));

mat_REF_6283, mat_SIF_6283, wlf_6283, wl_6283, df_6283 = derive_spectrum(df_6283, "LUT", 24)
mat_REF_6287, mat_SIF_6287, wlf_6287, wl_6287, df_6287 = derive_spectrum(df_6287, "LUT", 24)
mat_REF_6348, mat_SIF_6348, wlf_6348, wl_6348, df_6348 = derive_spectrum(df_6348, "LUT", 25)

site_name = "Niwot Ridge 2020-06-16 (CI, Cab, LAI)"

sif_plot = scatter_plot_SIF(df_6348, site_name)
annotate!(1.7, 0.65, text("R² = $(round(r2(reg740), digits = 2)) \nMAE = $mae_sif740 \n$pval740", :left, 10))
annotate!(1, 0.65, text("R² = $(round(r2(reg757), digits = 2)) \nMAE = $mae_sif757 \n$pval757", :left, 10))
annotate!(1, 0.1, text("R² = $(round(r2(reg771), digits = 2)) \nMAE = $mae_sif771 \n$pval771", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF_OCO3_CliMA_Scatter_Niwot_2020-06-16_CI.pdf")

ref_plot = scatter_plot_REF(df_6348, site_name)
annotate!(141, 95, text("R² = $(round(r2(reg757), digits = 2)) \nMAE = $mae_sif757 \n$pval757", :left, 10))
annotate!(133.5, 95, text("R² = $(round(r2(reg771), digits = 2)) \nMAE = $mae_sif771 \n$pval771", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/REF_OCO3_CliMA_Scatter_Niwot_2020-06-16_CI.pdf")