using CSV
code_dir = "C:/Russell/Projects/Geometry/Julia_Scripts/Code/"
include(code_dir*"Derive/derive_spectrum.jl");
include(code_dir*"Plot/Plot_OCO3_Data.jl");

df_7136 = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_umb_2020-08-06_7136.csv"
df_7136 = DataFrame!(CSV.File(df_7136));
df_7215 = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_umb_2020-08-11_7215.csv"
df_7215 = DataFrame!(CSV.File(df_7215));

# Filter data for ENF
# df_7136 = df_7136[df_7136[:IGBP_index] .== 1, :]
# df_7215 = df_7215[df_7215[:IGBP_index] .== 1, :]

# Filter data for ENF and Mixed Forest
# df_7136 = df_7136[findall(in([1,5]), df_7136[:IGBP_index]), :]
# df_7215 = df_7215[findall(in([1,5]), df_7215[:IGBP_index]), :]

mat_REF_7136_clump, mat_SIF_7136_clump, wlf_7136_clump, wl_7136_clump, df_7136_clump = derive_spectrum(df_7136, "LUT", 32, "clump", 1);
mat_REF_7215_clump, mat_SIF_7215_clump, wlf_7215_clump, wl_7215_clump, df_7215_clump = derive_spectrum(df_7215, "LUT", 33, "clump", 1);

mat_REF_7136_noclump, mat_SIF_7136_noclump, wlf_7136_noclump, wl_7136_noclump, df_7136_noclump = derive_spectrum(df_7136, "LUT", 32, 1, 1);
mat_REF_7215_noclump, mat_SIF_7215_noclump, wlf_7215_noclump, wl_7215_noclump, df_7215_noclump = derive_spectrum(df_7215, "LUT", 33, 1, 1);

site_name = "UMB August 6 and 11, 2020"

sif_results_clump, dfs_mean_clump = results_SIF([df_7136_clump, df_7215_clump], site_name)
sif_results_noclump, dfs_mean_noclump = results_SIF([df_7136_noclump, df_7215_noclump], site_name)

scatter(dfs_mean_clump.mean_sif757_sim, dfs_mean_clump.mean_sif757, xlabel = "Mean CliMA SIF₇₅₇ (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 SIF₇₅₇ (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
                title = site_name, titlefontsize = 13, legend = :topleft, label = "Clumping", framestyle = :box, yerror = dfs_mean_clump.stderr_sif757)
scatter!(dfs_mean_noclump.mean_sif757_sim, dfs_mean_noclump.mean_sif757, label = "No Clumping", reg = true, linewidth = 2, yerror = dfs_mean_noclump.stderr_sif757)
xlims!(0,1.5)
ylims!(0,1.5)
annotate!(0.6, 1.25, text(string("R² = ", sif_results_clump[1,2], "\nMAE = ", sif_results_clump[5,2], "\nslope = ", sif_results_clump[3,2], "\n", sif_results_clump[2,2]), :left, 10))
annotate!(1.1, 0.25, text(string("R² = ", sif_results_noclump[1,2], "\nMAE = ", sif_results_noclump[5,2], "\nslope = ", sif_results_noclump[3,2], "\n", sif_results_noclump[2,2]), :left, 10))

savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Renato/SIF_OCO3_CliMA_Scatter_UMB_Clumping.pdf")



ref_results, dfs_mean, ref_plot = scatter_plot_REF([df_7136, df_7215, df_6348], site_name)
annotate!(100, 45, text(string("R² = ", ref_results[1,1], "\nMAE = ", ref_results[5,1], "\n", ref_results[2,1]), :left, 10))
annotate!(70, 60, text(string("R² = ", ref_results[1,2], "\nMAE = ", ref_results[5,2], "\n", ref_results[2,2]), :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/REF_OCO3_CliMA_Scatter_Niwot_W_Clumping_Test.pdf")

CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_umb_means_noclump.csv", dfs_mean_noclump)
CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_umb_means_clump.csv", dfs_mean_clump)
# CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_2020-06-12_7136.csv", df_7136)
# CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_2020-06-12_7215.csv", df_7215)
# CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_2020-06-16_6348.csv", df_6348)

scatter(dfs_mean.mean_sif757, dfs_mean.mean_sif771, xlabel = "OCO-3 SIF 757nm (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "OCO-3 SIF 771nm (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
title = "Mean OCO-3 SIF Grouped by PA at Niwot 2020-06-12", titlefontsize = 13, legend = false)

scatter(dfs_mean.mean_sif757_sim, dfs_mean.mean_sif771_sim, xlabel = "CliMA SIF 757nm (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "CliMA SIF 771nm (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
title = "Mean CliMA SIF Grouped by PA at Niwot 2020-06-12", titlefontsize = 13, legend = false)