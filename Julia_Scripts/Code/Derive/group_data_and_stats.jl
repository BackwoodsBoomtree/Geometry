using Statistics
using GLM
using PlotPlants

# Inputs are vectors of dataframes
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
function stats_and_errors(dfs::Vector{DataFrame})
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

    # Build empty df for output of statistics
    r2_names        = ["r2_sif740","r2_sif757","r2_sif771","r2_ref757","r2_ref771"];
    pval_names      = ["pval_sif740","pval_sif757","pval_sif771","pval_ref757","pval_ref771"];
    slope_names     = ["slope_sif740","slope_sif757","slope_sif771","slope_ref757","slope_ref771"];
    intercept_names = ["intercept_sif740","intercept_sif757","intercept_sif771","intercept_ref757","intercept_ref771"];
    mae_names       = ["mae_sif740","mae_sif757","mae_sif771","mae_ref757","mae_ref771"];
    all_names       = vcat(r2_names, pval_names, slope_names, intercept_names, mae_names)
    stats_return    = DataFrame(fill(Float64, length(all_names)), Symbol.(all_names), 1);

    # SIF stats
    df740_sif, df757_sif, df771_sif    = DataFrame(X = mean_merge.mean_sif740_sim, Y = mean_merge.mean_sif740),
                                         DataFrame(X = mean_merge.mean_sif757_sim, Y = mean_merge.mean_sif757),
                                         DataFrame(X = mean_merge.mean_sif771_sim, Y = mean_merge.mean_sif771);
    reg740_sif, reg757_sif, reg771_sif = lm(@formula(Y ~ X), df740_sif), lm(@formula(Y ~ X), df757_sif), lm(@formula(Y ~ X), df771_sif);
    # Reflectance stats
    df757_ref, df771_ref               = DataFrame(X = mean_merge.mean_ref757_sim, Y = mean_merge.mean_ref757),
                                         DataFrame(X = mean_merge.mean_ref771_sim, Y = mean_merge.mean_ref771);
    reg757_ref, reg771_ref             = lm(@formula(Y ~ X), df757_ref), lm(@formula(Y ~ X), df771_ref);

    # Fill output df with SIF Stats
    stats_return.r2_sif740, stats_return.r2_sif757, stats_return.r2_sif771                      = # R2
        (round(r2(reg740_sif), digits = 2), round(r2(reg757_sif), digits = 2), round(r2(reg771_sif), digits = 2)) 
    stats_return.pval_sif740, stats_return.pval_sif757, stats_return.pval_sif771                = # Pval
        (round(coeftable(reg740_sif).cols[4][2], digits = 3), round(coeftable(reg757_sif).cols[4][2], digits = 3), round(coeftable(reg771_sif).cols[4][2], digits = 3)) 
    stats_return.slope_sif740, stats_return.slope_sif757, stats_return.slope_sif771             = # slope
        (round(coef(reg740_sif)[2], digits = 2), round(coef(reg757_sif)[2], digits = 2), round(coef(reg771_sif)[2], digits = 2))
    stats_return.intercept_sif740, stats_return.intercept_sif757, stats_return.intercept_sif771 = # intercept
        (round(coef(reg740_sif)[1], digits = 2), round(coef(reg757_sif)[1], digits = 2), round(coef(reg771_sif)[1], digits = 2))
    # Fill output df with REF Stats
    stats_return.r2_ref757, stats_return.r2_ref771               = # R2
        (round(r2(reg757_ref), digits = 2), round(r2(reg771_ref), digits = 2)) 
    stats_return.pval_ref757, stats_return.pval_ref771           = # Pval
        (round(coeftable(reg757_ref).cols[4][2], digits = 3), round(coeftable(reg771_ref).cols[4][2], digits = 3)) 
    stats_return.slope_ref757, stats_return.slope_ref771         = # slope
        (round(coef(reg757_ref)[2], digits = 2), round(coef(reg771_ref)[2], digits = 2))
    stats_return.intercept_ref757, stats_return.intercept_ref771 = # intercept
        (round(coef(reg757_ref)[1], digits = 2), round(coef(reg771_ref)[1], digits = 2))
    println("Computed regression stats for the ", nrow(mean_merge), " groups.")

    # MAE
    stats_return.mae_sif740 = round(PlotPlants.mae(mean_merge.mean_sif740, mean_merge.mean_sif740_sim), digits = 2)
    stats_return.mae_sif757 = round(PlotPlants.mae(mean_merge.mean_sif757, mean_merge.mean_sif757_sim), digits = 2)
    stats_return.mae_sif771 = round(PlotPlants.mae(mean_merge.mean_sif771, mean_merge.mean_sif771_sim), digits = 2)
    stats_return.mae_ref757 = round(PlotPlants.mae(mean_merge.mean_ref757, mean_merge.mean_ref757_sim), digits = 2)
    stats_return.mae_ref771 = round(PlotPlants.mae(mean_merge.mean_ref771, mean_merge.mean_ref771_sim), digits = 2)
    println("Calculated MAE for the ", nrow(mean_merge), " groups.")

    return(mean_merge, stats_return);
end