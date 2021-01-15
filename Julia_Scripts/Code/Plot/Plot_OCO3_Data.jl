code_dir = "C:/Russell/Projects/Geometry/Julia_Scripts/Code/"
include(code_dir*"Derive/group_data_and_stats.jl");
using Plots
using Plots.PlotMeasures

# Plotting SIF
function scatter_plot_SIF(dfs::Vector{DataFrame}, plot_title)
    dfs_mean, dfs_mae = means_and_errors(dfs)
    
    # df for regression results. Row 1 to 5 is R2, pval, slope, intercept, mae.
    sif_names = ["sif_740", "sif_757", "sif_771"]
    sif_results = DataFrame(fill(Any, length(sif_names)), Symbol.(sif_names), 5);
    
    # Scatterplot of group means for SIF 740, 757, 771
    p = scatter(dfs_mean.mean_sif740_sim, dfs_mean.mean_sif740, xlabel = "Mean CliMA SIF (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 SIF (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
                title = plot_title, titlefontsize = 13, legend = :topleft, label = "740 nm", framestyle = :box, yerror = dfs_mean.stderr_sif740)
    scatter!(dfs_mean.mean_sif757_sim, dfs_mean.mean_sif757, label = "757 nm", reg = true, linewidth = 2, yerror = dfs_mean.stderr_sif757)
    scatter!(dfs_mean.mean_sif771_sim, dfs_mean.mean_sif771, label = "771 nm", reg = true, linewidth = 2, yerror = dfs_mean.stderr_sif771)
    df740, df757, df771 = DataFrame(X = dfs_mean.mean_sif740_sim, Y = dfs_mean.mean_sif740),
                          DataFrame(X = dfs_mean.mean_sif757_sim, Y = dfs_mean.mean_sif757),
                          DataFrame(X = dfs_mean.mean_sif771_sim, Y = dfs_mean.mean_sif771)
    reg740, reg757, reg771 = lm(@formula(Y ~ X), df740), lm(@formula(Y ~ X), df757), lm(@formula(Y ~ X), df771)
    sif_results[1,:] = (round(r2(reg740), digits = 2), round(r2(reg757), digits = 2), round(r2(reg771), digits = 2)) # R2
    sif_results[2,:] = (round(coeftable(reg740).cols[4][2], digits = 3), round(coeftable(reg757).cols[4][2], digits = 3), round(coeftable(reg771).cols[4][2], digits = 3)) # Pval
    sif_results[3,:] = (round(coef(reg740)[2], digits = 2), round(coef(reg757)[2], digits = 2), round(coef(reg771)[2], digits = 2)) # slope
    sif_results[4,:] = (round(coef(reg740)[1], digits = 2), round(coef(reg757)[1], digits = 2), round(coef(reg771)[1], digits = 2)) # intercept
    sif_results[5,:] = (dfs_mae[1,1], dfs_mae[1,2], dfs_mae[1,3]) #mae
    for i in 1:3
        if sif_results[2,i] < 0.05
            sif_results[2,i] = "p-value ≤ 0.05"
        elseif sif_results[2,i] >= 0.05
            sif_results[2,i] = string("p-value = ", sif_results[2,i])
        end
    end
    return(sif_results, dfs_mean, p)
end

function scatter_plot_SIF_757_niwot(dfs::Vector{DataFrame}, plot_title)
    dfs_mean, dfs_mae = means_and_errors(dfs)
    
    # df for regression results. Row 1 to 5 is R2, pval, slope, intercept, mae.
    sif_names = ["sif_740", "sif_757", "sif_771"]
    sif_results = DataFrame(fill(Any, length(sif_names)), Symbol.(sif_names), 5);

    # Stats
    df740, df757, df771 = DataFrame(X = dfs_mean.mean_sif740_sim, Y = dfs_mean.mean_sif740),
                          DataFrame(X = dfs_mean.mean_sif757_sim, Y = dfs_mean.mean_sif757),
                          DataFrame(X = dfs_mean.mean_sif771_sim, Y = dfs_mean.mean_sif771)
    reg740, reg757, reg771 = lm(@formula(Y ~ X), df740), lm(@formula(Y ~ X), df757), lm(@formula(Y ~ X), df771)
    sif_results[1,:] = (round(r2(reg740), digits = 2), round(r2(reg757), digits = 2), round(r2(reg771), digits = 2)) # R2
    sif_results[2,:] = (round(coeftable(reg740).cols[4][2], digits = 3), round(coeftable(reg757).cols[4][2], digits = 3), round(coeftable(reg771).cols[4][2], digits = 3)) # Pval
    sif_results[3,:] = (round(coef(reg740)[2], digits = 2), round(coef(reg757)[2], digits = 2), round(coef(reg771)[2], digits = 2)) # slope
    sif_results[4,:] = (round(coef(reg740)[1], digits = 2), round(coef(reg757)[1], digits = 2), round(coef(reg771)[1], digits = 2)) # intercept
    sif_results[5,:] = (dfs_mae[1,1], dfs_mae[1,2], dfs_mae[1,3]) #mae
    for i in 1:3
        if sif_results[2,i] < 0.05
            sif_results[2,i] = "p-value ≤ 0.05"
        elseif sif_results[2,i] >= 0.05
            sif_results[2,i] = string("p-value = ", sif_results[2,i])
        end
    end
    return(sif_results, dfs_mean)
end

# Plotting SIF relative
function scatter_plot_SIF_relative(dfs::Vector{DataFrame}, plot_title)
    dfs_mean, dfs_mae = means_and_errors(dfs)
    
    # df for regression results. Row 1 to 5 is R2, pval, slope, intercept, mae.
    sif_names = ["sif_740", "sif_757", "sif_771"]
    sif_results = DataFrame(fill(Any, length(sif_names)), Symbol.(sif_names), 5);
    
    # Scatterplot of group means for SIF 740, 757, 771
    p = scatter(dfs_mean.mean_sif740_sim, dfs_mean.mean_sif740, xlabel = "Mean CliMA SIF (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 SIF (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
                title = plot_title, titlefontsize = 13, legend = :topleft, label = "740 nm", framestyle = :box, yerror = dfs_mean.stderr_sif740)
    scatter!(dfs_mean.mean_sif757_sim, dfs_mean.mean_sif757, label = "757 nm", reg = true, linewidth = 2, yerror = dfs_mean.stderr_sif757)
    scatter!(dfs_mean.mean_sif771_sim, dfs_mean.mean_sif771, label = "771 nm", reg = true, linewidth = 2, yerror = dfs_mean.stderr_sif771)
    df740, df757, df771 = DataFrame(X = dfs_mean.mean_sif740_sim, Y = dfs_mean.mean_sif740),
                          DataFrame(X = dfs_mean.mean_sif757_sim, Y = dfs_mean.mean_sif757),
                          DataFrame(X = dfs_mean.mean_sif771_sim, Y = dfs_mean.mean_sif771)
    reg740, reg757, reg771 = lm(@formula(Y ~ X), df740), lm(@formula(Y ~ X), df757), lm(@formula(Y ~ X), df771)
    sif_results[1,:] = (round(r2(reg740), digits = 2), round(r2(reg757), digits = 2), round(r2(reg771), digits = 2)) # R2
    sif_results[2,:] = (round(coeftable(reg740).cols[4][2], digits = 3), round(coeftable(reg757).cols[4][2], digits = 3), round(coeftable(reg771).cols[4][2], digits = 3)) # Pval
    sif_results[3,:] = (round(coef(reg740)[2], digits = 2), round(coef(reg757)[2], digits = 2), round(coef(reg771)[2], digits = 2)) # slope
    sif_results[4,:] = (round(coef(reg740)[1], digits = 2), round(coef(reg757)[1], digits = 2), round(coef(reg771)[1], digits = 2)) # intercept
    sif_results[5,:] = (dfs_mae[1,1], dfs_mae[1,2], dfs_mae[1,3]) #mae
    for i in 1:3
        if sif_results[2,i] < 0.05
            sif_results[2,i] = "p-value ≤ 0.05"
        elseif sif_results[2,i] >= 0.05
            sif_results[2,i] = string("p-value = ", sif_results[2,i])
        end
    end
    return(sif_results, dfs_mean, p)
end

# Plotting REF
function scatter_plot_REF(dfs::Vector{DataFrame}, plot_title)
    dfs_mean, dfs_mae = means_and_errors(dfs)

    # df for regression results. Row 1 to 5 is R2, pval, slope, intercept, mae.
    ref_names = ["ref_757", "ref_771"]
    ref_results = DataFrame(fill(Any, length(ref_names)), Symbol.(ref_names), 5);

    # Scatterplot of group means for reflected radiance at 757 and 771
    p = scatter(dfs_mean.mean_ref757_sim, dfs_mean.mean_ref757, xlabel = "Mean CliMA Reflected Radiance (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 Reflected Radiance (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
                title = plot_title, titlefontsize = 13, legend = :topleft, label = "757 nm", framestyle = :box, yerror = dfs_mean.stderr_ref757)
    scatter!(dfs_mean.mean_ref771_sim, dfs_mean.mean_ref771, label = "771 nm", reg = true, linewidth = 2, yerror = dfs_mean.stderr_ref771)
    df757, df771 = DataFrame(X = dfs_mean.mean_ref757_sim, Y = dfs_mean.mean_ref757), DataFrame(X = dfs_mean.mean_ref771_sim, Y = dfs_mean.mean_ref771)
    reg757, reg771 = lm(@formula(Y ~ X), df757), lm(@formula(Y ~ X), df771)
    ref_results[1,:] = (round(r2(reg757), digits = 2), round(r2(reg771), digits = 2)) # R2
    ref_results[2,:] = (round(coeftable(reg757).cols[4][2], digits = 3), round(coeftable(reg771).cols[4][2], digits = 3)) # Pval
    ref_results[3,:] = (round(coef(reg757)[2], digits = 2), round(coef(reg771)[2], digits = 2)) # slope
    ref_results[4,:] = (round(coef(reg757)[1], digits = 2), round(coef(reg771)[1], digits = 2)) # intercept
    ref_results[5,:] = (dfs_mae[1,4], dfs_mae[1,5]) #mae
    for i in 1:2
        if ref_results[2,i] < 0.05
            ref_results[2,i] = "p-value ≤ 0.05"
        elseif ref_results[2,i] >= 0.05
            ref_results[2,i] = string("p-value = ", ref_results[2,i])
        end
    end
    return(ref_results, dfs_mean, p)
end

