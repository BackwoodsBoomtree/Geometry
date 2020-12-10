using PlotPlants
using Statistics
using GLM
using Plots
using Plots.PlotMeasures

# Determines groups of data by phase angle. Adds a group_by_pa column and populates with group number
function group_by_pa(dfs::Vector{DataFrame})
    dfs_return = Vector{DataFrame}()
    for k in 1:length(dfs)
        if "group_by_pa" ∉ names(dfs[k])
            insertcols!(dfs[k], ncol(dfs[k])+1, :group_by_pa => NaN)
        end
        sort!(dfs[k], :PA)
        global j = 1
        for i in 1:nrow(dfs[k])
            if i == 1
                dfs[k].group_by_pa[i] = j
            elseif i != 1
                if abs(dfs[k].PA[i - 1] - dfs[k].PA[i]) < 0.25
                    dfs[k].group_by_pa[i] = j
                else
                    global j = j + 1
                    dfs[k].group_by_pa[i] = j
                end
            end
        end
        push!(dfs_return, dfs[k])
    end
    return(dfs_return)
end

# For each dataframe, calculates means and standard errors for each group, and returns dataframe of such for each input dataframe
function means_and_errors(dfs::Vector{DataFrame})
    group_by_pa(dfs)
    mean_return = Vector{DataFrame}()
    mean_names = ["mean_sif740", "mean_sif740_sim", "stderr_sif740",
                 "mean_sif757", "mean_sif757_sim", "stderr_sif757",
                 "mean_sif771", "mean_sif771_sim", "stderr_sif771",
                 "mean_ref757", "mean_ref757_sim", "stderr_ref757",
                 "mean_ref771", "mean_ref771_sim", "stderr_ref771",
                 "mean_raa", "mean_vza"]
    for k in 1:length(dfs)
        # Build temp df with rows equal to number of groups, but temp df rows that don't get populated with minimum number
        # n in groups gets removed
        println("Dataframe ", k)
        dfs_temp = DataFrame(fill(Float64, length(mean_names)), Symbol.(mean_names), length(unique(dfs[k].group_by_pa)))
        for col in eachcol(dfs_temp)
            col .= NaN
        end
        for i in 1:length(unique(dfs[k].group_by_pa))
            if length(dfs[k][dfs[k][:group_by_pa] .== [i], :].group_by_pa) > 9
                dfs_temp.mean_sif740[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].SIF_740nm)
                dfs_temp.mean_sif740_sim[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].SIF740_Sim)
                dfs_temp.stderr_sif740[i] = std(dfs[k][dfs[k][:group_by_pa] .== [i], :].SIF_740nm) / sqrt(nrow(dfs[k][dfs[k][:group_by_pa] .== [i], :]))

                dfs_temp.mean_sif757[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].SIF_757nm)
                dfs_temp.mean_sif757_sim[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].SIF757_Sim)
                dfs_temp.stderr_sif757[i] = std(dfs[k][dfs[k][:group_by_pa] .== [i], :].SIF_757nm) / sqrt(nrow(dfs[k][dfs[k][:group_by_pa] .== [i], :]))

                dfs_temp.mean_sif771[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].SIF_771nm)
                dfs_temp.mean_sif771_sim[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].SIF771_Sim)
                dfs_temp.stderr_sif771[i] = std(dfs[k][dfs[k][:group_by_pa] .== [i], :].SIF_771nm) / sqrt(nrow(dfs[k][dfs[k][:group_by_pa] .== [i], :]))

                dfs_temp.mean_ref757[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].continuum_radiance_757nm)
                dfs_temp.mean_ref757_sim[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].REF757_Sim)
                dfs_temp.stderr_ref757[i] = std(dfs[k][dfs[k][:group_by_pa] .== [i], :].continuum_radiance_757nm) / sqrt(nrow(dfs[k][dfs[k][:group_by_pa] .== [i], :]))

                dfs_temp.mean_ref771[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].continuum_radiance_771nm)
                dfs_temp.mean_ref771_sim[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].REF771_Sim)
                dfs_temp.stderr_ref771[i] = std(dfs[k][dfs[k][:group_by_pa] .== [i], :].continuum_radiance_771nm) / sqrt(nrow(dfs[k][dfs[k][:group_by_pa] .== [i], :]))

                dfs_temp.mean_raa[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].RAz)
                dfs_temp.mean_vza[i] = mean(dfs[k][dfs[k][:group_by_pa] .== [i], :].VZA)

                println("Group ", i, " included and pushed!")
            else
                println("Group ", i, " excluded with an N of ", length(dfs[k][dfs[k][:group_by_pa] .== [i], :].group_by_pa), " from dataframe ", k)
            end
        end
        dfs_temp = filter(:mean_sif740 => mean_sif740 -> !(isnan(mean_sif740)), dfs_temp) # filter out rows that didn't get populated
        push!(mean_return, dfs_temp)
        println("Number of groups included: ", length(dfs_temp.mean_sif740), " for dataframe ", k)
    end
    # Append dfs from mean_return for output so we can calculate mae
    if length(mean_return) > 1
        for df in 1:length(mean_return) - 1
            if df == 1
                mean_merge = vcat(mean_return[df], mean_return[df+1]);
            else
                append!(mean_merge, mean_return[df+1]);
            end
        end
    else
        mean_merge = mean_return[1];
    end
    mean_return = nothing # Don't need it any more
    println("Calculated MAE for the ", nrow(mean_merge), " groups.")
    mae_names = ["mae_sif740","mae_sif757","mae_sif771","mae_ref757","mae_ref771"]
    mae_return = DataFrame(fill(Float64, length(mae_names)), Symbol.(mae_names), 1);
    mae_return.mae_sif740 = round(PlotPlants.mae(mean_merge.mean_sif740, mean_merge.mean_sif740_sim), digits = 2)
    mae_return.mae_sif757 = round(PlotPlants.mae(mean_merge.mean_sif757, mean_merge.mean_sif757_sim), digits = 2)
    mae_return.mae_sif771 = round(PlotPlants.mae(mean_merge.mean_sif771, mean_merge.mean_sif771_sim), digits = 2)
    mae_return.mae_ref757 = round(PlotPlants.mae(mean_merge.mean_ref757, mean_merge.mean_ref757_sim), digits = 2)
    mae_return.mae_ref771 = round(PlotPlants.mae(mean_merge.mean_ref771, mean_merge.mean_ref771_sim), digits = 2)

    return(mean_merge, mae_return);
end

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

