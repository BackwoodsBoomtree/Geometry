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

# Filter data for landcover
# df_6283 = df_6283[df_6283[:IGBP_index] .== 1, :]
# df_6287 = df_6287[df_6287[:IGBP_index] .== 1, :]
# df_6348 = df_6348[df_6348[:IGBP_index] .== 1, :]


mat_REF_6283, mat_SIF_6283, wlf_6283, wl_6283, df_6283 = derive_spectrum(df_6283, "LUT", 24);
mat_REF_6287, mat_SIF_6287, wlf_6287, wl_6287, df_6287 = derive_spectrum(df_6287, "LUT", 24);
mat_REF_6348, mat_SIF_6348, wlf_6348, wl_6348, df_6348 = derive_spectrum(df_6348, "LUT", 25);

site_name = "Niwot Ridge June 12 and 16 (with clumping index)"

sif_results, dfs_mean, sif_plot = scatter_plot_SIF([df_6283, df_6287, df_6348], site_name)
annotate!(1.75, 0.2, text(string("R² = ", sif_results[1,1], "\nMAE = ", sif_results[5,1], "\n", sif_results[2,1]), :left, 10))
annotate!(1, 0.65, text(string("R² = ", sif_results[1,2], "\nMAE = ", sif_results[5,2], "\n", sif_results[2,2]), :left, 10))
annotate!(0.4, 0.35, text(string("R² = ", sif_results[1,3], "\nMAE = ", sif_results[5,3], "\n", sif_results[2,3]), :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF_OCO3_CliMA_Scatter_Niwot_W_Clumping_Test.pdf")

ref_results, dfs_mean, ref_plot = scatter_plot_REF([df_6283, df_6287, df_6348], site_name)
annotate!(100, 45, text(string("R² = ", ref_results[1,1], "\nMAE = ", ref_results[5,1], "\n", ref_results[2,1]), :left, 10))
annotate!(70, 60, text(string("R² = ", ref_results[1,2], "\nMAE = ", ref_results[5,2], "\n", ref_results[2,2]), :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/REF_OCO3_CliMA_Scatter_Niwot_W_Clumping_Test.pdf")

CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_means_noclump.csv", dfs_mean)
# CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_2020-06-12_6283.csv", df_6283)
# CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_2020-06-12_6287.csv", df_6287)
# CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_2020-06-16_6348.csv", df_6348)

scatter(dfs_mean.mean_sif757, dfs_mean.mean_sif771, xlabel = "OCO-3 SIF 757nm (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "OCO-3 SIF 771nm (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
title = "Mean OCO-3 SIF Grouped by PA at Niwot 2020-06-12", titlefontsize = 13, legend = false)

scatter(dfs_mean.mean_sif757_sim, dfs_mean.mean_sif771_sim, xlabel = "CliMA SIF 757nm (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "CliMA SIF 771nm (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
title = "Mean CliMA SIF Grouped by PA at Niwot 2020-06-12", titlefontsize = 13, legend = false)