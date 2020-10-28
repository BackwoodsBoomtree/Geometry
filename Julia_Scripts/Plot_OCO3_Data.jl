using PlotPlants
using Statistics
using GLM
using Plots
using Plots.PlotMeasures


df1 = DataFrame(REF = 1:250, PA=rand(250))
df2 = DataFrame(REF = 1:250, PA=rand(250))
df3 = DataFrame(REF = 1:250, PA=rand(250))
dfs_new = group_by_pa([df1, df2, df3])

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
    dfs_return = Vector{DataFrame}()
    col_names = ["mean_sif740", "mean_sif740_sim", "stderr_sif740",
                 "mean_sif757", "mean_sif757_sim", "stderr_sif757",
                 "mean_sif771", "mean_sif771_sim", "stderr_sif771",
                 "mean_ref757", "mean_ref757_sim", "stderr_ref757",
                 "mean_ref771", "mean_ref771_sim", "stderr_ref771",
                 "mean_raa", "mean_vza"]
    # # Means and Standard Errors for each group_by_pa
    # global mean_sif740, mean_sif740_sim, stderr_sif740  = Float64[], Float64[], Float64[]
    # global mean_sif757, mean_sif757_sim, stderr_sif757  = Float64[], Float64[], Float64[]
    # global mean_sif771, mean_sif771_sim, stderr_sif771  = Float64[], Float64[], Float64[]
    # global mean_ref757, mean_ref757_sim, stderr_ref757  = Float64[], Float64[], Float64[]
    # global mean_ref771, mean_ref771_sim, stderr_ref771  = Float64[], Float64[], Float64[]
    # global mean_raa, mean_vza  = Float64[], Float64[]
    for k in 1:length(dfs)
        # Build temp df with rows equal to number of groups, but temp df rows that don't get populated with minimum number
        # n in groups gets removed
        println("Dataframe ", k)
        dfs_temp = DataFrame(fill(Float64, length(col_names)), Symbol.(col_names), length(unique(dfs[k].group_by_pa)))
        println(dfs_temp)
        for i in 1:length(unique(dfs[k].group_by_pa))
            println("For Group ", i)
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

                println(dfs_temp)
                push!(dfs_return, dfs_temp)
                println("Group ", i, " included and pushed!")
            else
                println("Group ", i, " excluded with an N of ", length(dfs[k][dfs[k][:group_by_pa] .== [i], :].group_by_pa), " from dataframe ", k)
            end
        end
        println("Number of groups includedadfasfwreqr4343: ", length(mean_sif740), " for dataframe ", k)
    end
    dfs_temp <- filter(:mean_sif740 => mean_sif740 -> !(isnan(mean_sif740)), dfs_temp) # filter out rows that didn't get populated
    return(dfs_return)
end

dfs_test = means_and_errors([df_6283])

# dfs_new_means = means_and_errors(dfs_new)



# # Mean aboslute error
# global mae_sif740 = round(PlotPlants.mae(mean_sif740, mean_sif740_sim), digits = 2)
# global mae_sif757 = round(PlotPlants.mae(mean_sif757, mean_sif757_sim), digits = 2)
# global mae_sif771 = round(PlotPlants.mae(mean_sif771, mean_sif771_sim), digits = 2)
# global mae_ref757 = round(PlotPlants.mae(mean_ref757, mean_ref757_sim), digits = 2)
# global mae_ref771 = round(PlotPlants.mae(mean_ref771, mean_ref771_sim), digits = 2)

function scatter_plot_SIF(dfs::Vector{DataFrame}, plot_title)
    group_by_pa(df)
    means_and_errors(df)
    
    # Scatterplot of group means for SIF 740, 757, 771
    p = scatter(mean_sif740_sim, mean_sif740, xlabel = "Mean CliMA SIF (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 SIF (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
                title = plot_title, titlefontsize = 13, legend = :topleft, label = "740 nm", framestyle = :box, yerror = stderr_sif740)
    scatter!(mean_sif757_sim, mean_sif757, label = "757 nm", reg = true, linewidth = 2, yerror = stderr_sif757)
    scatter!(mean_sif771_sim, mean_sif771, label = "771 nm", reg = true, linewidth = 2, yerror = stderr_sif771)
    df740, df757, df771 = DataFrame(X = mean_sif740_sim, Y = mean_sif740), DataFrame(X = mean_sif757_sim, Y = mean_sif757), DataFrame(X = mean_sif771_sim, Y = mean_sif771)
    global reg740, reg757, reg771 = lm(@formula(Y ~ X), df740), lm(@formula(Y ~ X), df757), lm(@formula(Y ~ X), df771)
    slope740, slope757, slope717 = round(coef(reg740)[2], digits = 2), round(coef(reg757)[2], digits = 2), round(coef(reg771)[2], digits = 2)
    intercept740, intercept757, intercept771 = round(coef(reg740)[1]), round(coef(reg757)[1]), round(coef(reg771)[1])
    pval740, pval757, pval771 = round(coeftable(reg740).cols[4][2], digits = 3), round(coeftable(reg757).cols[4][2], digits = 3), round(coeftable(reg771).cols[4][2], digits = 3)
    if pval757 == 0
        global pval757 = "p-value ≤ 0.001"
        global pval771 = "p-value ≤ 0.001"
        global pval740 = "p-value ≤ 0.001"
    else
        global pval757 = "p-value = $pval757"
        global pval771 = "p-value = $pval771"
        global pval740 = "p-value = $pval740"
    end
    return p
end

function scatter_plot_REF(df, plot_title)
    group_by_pa(df)
    means_and_errors(df)

    # Scatterplot of group means for reflected radiance at 757 and 771
    p = scatter(mean_ref757_sim, mean_ref757, xlabel = "Mean CliMA Reflected Radiance (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 Reflected Radiance (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
    title = plot_title, titlefontsize = 13, legend = :topleft, label = "757 nm", framestyle = :box, yerror = stderr_ref757)
    scatter!(mean_ref771_sim, mean_ref771, label = "771 nm", reg = true, linewidth = 2, yerror = stderr_ref771)
    df757, df771 = DataFrame(X = mean_ref757_sim, Y = mean_ref757), DataFrame(X = mean_ref771_sim, Y = mean_ref771)
    global reg757, reg771 = lm(@formula(Y ~ X), df757), lm(@formula(Y ~ X), df771)
    global slope757, slope717 = round(coef(reg757)[2], digits = 2), round(coef(reg771)[2], digits = 2)
    global intercept757, intercept771 = round(coef(reg757)[1]), round(coef(reg771)[1])
    pval757, pval771 = round(coeftable(reg757).cols[4][2], digits = 3), round(coeftable(reg771).cols[4][2], digits = 3)
    if pval757 == 0
        global pval757 = "p-value ≤ 0.001"
        global pval771 = "p-value ≤ 0.001"
    else
        global pval757 = "p-value = $pval757"
        global pval771 = "p-value = $pval771"
    end
    return p
end

