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

# Filter data for ENF
# df_6283 = df_6283[df_6283[:IGBP_index] .== 1, :]
# df_6287 = df_6287[df_6287[:IGBP_index] .== 1, :]
# df_6348 = df_6348[df_6348[:IGBP_index] .== 1, :]

# Filter data for ENF and Mixed Forest
# df_6283 = df_6283[findall(in([1,5]), df_6283[:IGBP_index]), :]
# df_6287 = df_6287[findall(in([1,5]), df_6287[:IGBP_index]), :]
# df_6348 = df_6348[findall(in([1,5]), df_6348[:IGBP_index]), :]

mat_REF_6283_clump, mat_SIF_6283_clump, wlf_6283_clump, wl_6283_clump, df_6283_clump = derive_spectrum(df_6283, "LUT", 24, "clump");
mat_REF_6287_clump, mat_SIF_6287_clump, wlf_6287_clump, wl_6287_clump, df_6287_clump = derive_spectrum(df_6287, "LUT", 24, "clump");
mat_REF_6348_clump, mat_SIF_6348_clump, wlf_6348_clump, wl_6348_clump, df_6348_clump = derive_spectrum(df_6348, "LUT", 25, "clump");

mat_REF_6283_noclump, mat_SIF_6283_noclump, wlf_6283_noclump, wl_6283_noclump, df_6283_noclump = derive_spectrum(df_6283, "LUT", 24, "noclump");
mat_REF_6287_noclump, mat_SIF_6287_noclump, wlf_6287_noclump, wl_6287_noclump, df_6287_noclump = derive_spectrum(df_6287, "LUT", 24, "noclump");
mat_REF_6348_noclump, mat_SIF_6348_noclump, wlf_6348_noclump, wl_6348_noclump, df_6348_noclump = derive_spectrum(df_6348, "LUT", 25, "noclump");

site_name = "Niwot Ridge June 12 and 16"

sif_results_clump, dfs_mean_clump = scatter_plot_SIF_757_niwot([df_6283_clump, df_6287_clump, df_6348_clump], site_name)
sif_results_noclump, dfs_mean_noclump = scatter_plot_SIF_757_niwot([df_6283_noclump, df_6287_noclump, df_6348_noclump], site_name)

scatter(dfs_mean_clump.mean_sif757_sim, dfs_mean_clump.mean_sif757, xlabel = "Mean CliMA SIF₇₅₇ (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 SIF₇₅₇ (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
                title = site_name, titlefontsize = 13, legend = :topleft, label = "Clumping", framestyle = :box, yerror = dfs_mean_clump.stderr_sif757)
scatter!(dfs_mean_noclump.mean_sif757_sim, dfs_mean_noclump.mean_sif757, label = "No Clumping", reg = true, linewidth = 2, yerror = dfs_mean_noclump.stderr_sif757)
xlims!(0.3,1)
ylims!(0.1,0.7)
annotate!(0.35, 0.525, text(string("R² = ", sif_results_clump[1,2], "\nMAE = ", sif_results_clump[5,2], "\nslope = ", sif_results_clump[3,2], "\n", sif_results_clump[2,2]), :left, 10))
annotate!(0.8, 0.2, text(string("R² = ", sif_results_noclump[1,2], "\nMAE = ", sif_results_noclump[5,2], "\nslope = ", sif_results_noclump[3,2], "\n", sif_results_noclump[2,2]), :left, 10))

savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF_OCO3_CliMA_Scatter_Niwot_W_Clumping_Test.pdf")

ref_results, dfs_mean, ref_plot = scatter_plot_REF([df_6283, df_6287, df_6348], site_name)
annotate!(100, 45, text(string("R² = ", ref_results[1,1], "\nMAE = ", ref_results[5,1], "\n", ref_results[2,1]), :left, 10))
annotate!(70, 60, text(string("R² = ", ref_results[1,2], "\nMAE = ", ref_results[5,2], "\n", ref_results[2,2]), :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/REF_OCO3_CliMA_Scatter_Niwot_W_Clumping_Test.pdf")

CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_means_noclump.csv", dfs_mean_noclump)
CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_means_clump.csv", dfs_mean_clump)
# CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_2020-06-12_6283.csv", df_6283)
# CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_2020-06-12_6287.csv", df_6287)
# CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_2020-06-16_6348.csv", df_6348)

scatter(dfs_mean.mean_sif757, dfs_mean.mean_sif771, xlabel = "OCO-3 SIF 757nm (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "OCO-3 SIF 771nm (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
title = "Mean OCO-3 SIF Grouped by PA at Niwot 2020-06-12", titlefontsize = 13, legend = false)

scatter(dfs_mean.mean_sif757_sim, dfs_mean.mean_sif771_sim, xlabel = "CliMA SIF 757nm (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "CliMA SIF 771nm (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
title = "Mean CliMA SIF Grouped by PA at Niwot 2020-06-12", titlefontsize = 13, legend = false)