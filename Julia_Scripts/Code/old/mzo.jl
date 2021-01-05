using CSV
include("Build_OCO3_Data.jl")
include("derive_spectrum.jl")
include("Plot_OCO3_Data.jl")

list_files = filter(x->contains(x, ".nc"), readdir("C:/Russell/Projects/Geometry/Data/oco3/Kessler", join = true))

df = build_oco3_data(list_files)

df = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_Ozark_USA_2020-08-20_7353.csv"
df = DataFrame!(CSV.File(df));

mat_REF, mat_SIF, wlf, wl, df = derive_spectrum(df)

sif_plot = scatter_plot_SIF(df, "Ozark USA")
annotate!(2.75, 2.75, text("R² = $(round(r2(reg740), digits = 2)) \nMAE = $mae_sif740 \n$pval740", :left, 10))
annotate!(2.6, 1.25, text("R² = $(round(r2(reg757), digits = 2)) \nMAE = $mae_sif757 \n$pval757", :left, 10))
annotate!(1.5, 1.0, text("R² = $(round(r2(reg771), digits = 2)) \nMAE = $mae_sif771 \n$pval771", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF_OCO3_CliMA_Scatter_Ozark_2020-08-20.pdf")

ref_plot = scatter_plot_REF(df, "Ozark USA")
annotate!(141, 95, text("R² = $(round(r2(reg757), digits = 2)) \nMAE = $mae_sif757 \n$pval757", :left, 10))
annotate!(133.5, 95, text("R² = $(round(r2(reg771), digits = 2)) \nMAE = $mae_sif771 \n$pval771", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/REF_OCO3_CliMA_Scatter_Ozark_2020-08-20.pdf")

# CSV.write("C:/Russell/Projects/Geometry/Data/model_output/Kessler.csv", df)