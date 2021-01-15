code_dir = "C:/Russell/Projects/Geometry/Julia_Scripts/Code/"
include(code_dir*"Derive/derive_spectrum.jl");
include(code_dir*"Derive/group_data_and_stats.jl");

using CSV
using Plots

amazon_data_path   = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_ATTO_Tower_Manaus_Brazil_(incorrect)_2020-06-26_6500.csv";
amazon_cab         = 60;
amazon_week        = NaN;
amazon_clumping    = 1;
amazon_sif_yield   = 0.5;
amazon_file_name   = "C:/Russell/Projects/Geometry/Julia_Scripts/CSV/ATTO_(incorrect)_2020-06-26_6500";

# Run model, get means and maes of grouped data
amazon_mat_REF, amazon_mat_SIF, wlf, wl, 
amazon_data   = derive_spectrum(amazon_data_path, amazon_cab, amazon_week, amazon_clumping, amazon_sif_yield);
amazon_mean, 
amazon_stats  = stats_and_errors([amazon_data]);

# Save model output and stats
CSV.write(amazon_file_name*"_simulated.csv", amazon_data);
CSV.write(amazon_file_name*"_means.csv", amazon_mean);
CSV.write(amazon_file_name*"_stats.csv", amazon_stats);

# Niwot Data
niwot_mean_clump = DataFrame!(CSV.File("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_means_clump.csv"));
niwot_mean_noclump = DataFrame!(CSV.File("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/OCO3_CliMA_niwot_means_noclump.csv"));

# Plot
scalefontsizes(1.3)
# scalefontsizes(1 / 1.3) # reset

main_layout = @layout([a b])

# Amazon
p_amazon = scatter(amazon_mean.mean_sif757_sim, amazon_mean.mean_sif757, xlabel = "Mean CliMA SIF (mW m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 SIF (mW m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
            title = "ATTO (incorrect)", titlefontsize = 13, legend = :topleft, label = "757 nm", framestyle = :box, yerror = amazon_mean.stderr_sif757)
scatter!(amazon_mean.mean_sif771_sim, amazon_mean.mean_sif771, label = "771 nm", reg = true, linewidth = 2, yerror = amazon_mean.stderr_sif771)
if amazon_stats.pval_sif757[1] < 0.05
    pval_sif757 = "p-value < 0.05"
elseif amazon_stats.pval_sif757[1] >= 0.05
    pval_sif757 = string("p-value = ", amazon_stats.pval_sif757)
end
if amazon_stats.pval_sif771[1] < 0.05
    pval_sif771 = "p-value < 0.05"
elseif amazon_stats.pval_sif771[1] >= 0.05
    pval_sif771 = string("p-value = ", amazon_stats.pval_sif771)
end
annotate!(1.15, 1.1, text(string("R² = ", amazon_stats.r2_sif757[1], "\nMAE = ", amazon_stats.mae_sif757[1], "\n", pval_sif757), :left, 12))
annotate!(0.80, 1.1, text(string("R² = ", amazon_stats.r2_sif771[1], "\nMAE = ", amazon_stats.mae_sif771[1], "\n", pval_sif771), :left, 12))

# Niwot
p_niwot = scatter(niwot_mean_clump.mean_sif757_sim, niwot_mean_clump.mean_sif757, xlabel = "Mean CliMA SIF₇₅₇ (mW m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 SIF₇₅₇ (mW m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
            title = "Niwot Ridge", titlefontsize = 13, legend = :topleft, label = "Clumping", framestyle = :box, yerror = niwot_mean_clump.stderr_sif757, color = palette(:default)[3])
scatter!(niwot_mean_noclump.mean_sif757_sim, niwot_mean_noclump.mean_sif757, label = "No Clumping", reg = true, linewidth = 2, yerror = niwot_mean_noclump.stderr_sif757, color = palette(:default)[4])
# Don't think I saved stats, but they were posted to Renato in slack
annotate!(0.4, 0.49, text(string("R² = 0.58", "\nMAE = 0.18", "\n", "p-value < 0.05"), :left, 12))
annotate!(0.8, 0.25, text(string("R² = 0.46", "\nMAE = 0.40", "\n", "p-value < 0.05"), :left, 12))

p_main = plot(p_amazon, p_niwot, layout = main_layout, size = (900,450))

savefig(p_main, "C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Proposal/Proposal_Plot_Amazon_Niwot.pdf")