using CSV
using DataFrames
using CanopyLayers
using PlotPlants
using Plots
using Plots.PlotMeasures
using GLM
using Statistics

struct SIFComparison2020{FT} end

"""
    derive_spectrum()
Derive and plot the albedo and SIF spectrum
"""
function derive_spectrum()
    FT = Float32;

    # read data from file
    # TODO add clumping factor
    # TODO add cab
    # replace the file name with your own file
    #_file = joinpath(@__DIR__, "../../data/russ_data.csv");
    _file = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_ecostress_us_syv_2020-06-17.csv"
    _data = DataFrame!(CSV.File(_file));

    # initialize canopy radiation module
    angles, can, can_opt, can_rad, in_rad,
    leaves, rt_con, rt_dim, soil, wls = initialize_rt_module(FT);

    # create a copy of incoming radiation
    in_rad_bak = deepcopy(in_rad);
    e_all_dire = sum(in_rad_bak.E_direct  .* wls.dWL) / 1000;
    e_all_diff = sum(in_rad_bak.E_diffuse .* wls.dWL) / 1000;

    # create a matrix to store the spectrum
    mat_REF1 = zeros(FT, (length(_data.vza), length(wls.WL)));
    mat_REF = zeros(FT, (length(_data.vza), length(wls.WL)));
    mat_SIF = zeros(FT, (length(_data.vza), length(wls.WLF)));

    # iterate through the data
    for i in eachindex(_data.vza)
        # change the angles
        angles.tts = _data.sza[i];
        angles.psi = _data.raa[i];
        angles.tto = _data.vza[i];

        # change the canopy profiles
        can.LAI    = _data.lai_cop[i];
        can.iLAI   = _data.lai_cop[i] * can.Ω / can.nLayer;

        # static LAI
        # can.LAI    = 5.7;
        # can.iLAI   = 5.7 * can.Ω / can.nLayer;

        # change incoming radiation
        in_rad = deepcopy(in_rad_bak);
        in_rad.E_direct  .*= _data.incoming_direct_era5[i]  / e_all_dire;
        in_rad.E_diffuse .*= _data.incoming_diffuse_era5[i] / e_all_diff;

        # Run fluspect on each layer to change Cab
        for i in 1:20
            leaves[i].Cab = 40;
            fluspect!(leaves[i], wls);
        end

        # re-run the simulations
        canopy_geometry!(can, angles, can_opt, rt_con);
        canopy_matrices!(leaves, can_opt);
        short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
        canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
        SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

        mat_REF[i,:] .= can_rad.Lo;
        mat_SIF[i,:] .= can_rad.SIF_obs;
    end

    # 742 nm and 737 nm to 740 nm
    # 757 nm to 757 nm
    # 767 nm and 774.5 nm to 771 nm
    _data.sif740_sim = mat_SIF[:,19] .* 0.4 .+ mat_SIF[:,20] .* 0.6;
    _data.sif757_sim = mat_SIF[:,23];
    _data.sif771_sim = mat_SIF[:,25] .* 0.4667 .+ mat_SIF[:,26] .* 0.5333;
    _data.ref757_sim = mat_REF[:,47];
    _data.ref771_sim = mat_REF[:,49] .* 0.4667 .+ mat_REF[:,50] .* 0.5333;

    return mat_REF, mat_SIF, wls.WLF, wls.WL, _data
    #return nothing
end

# Run model
mat_REF, mat_SIF, wlf, wl, data = derive_spectrum()

# Remove low LAI
data = data[(data[:lai_cop].>4), :]

# Remove high sza
data = data[(data[:sza].<32), :]

# Groups
if "group" in names(data) == false
    insertcols!(data, ncol(data)+1, :group => NaN)
end
sort!(data, :pa)
j = 1
for i in 1:nrow(data)
    if i == 1
        data.group[i] = j
    elseif i != 1
        if abs(data.pa[i - 1] - data.pa[i]) < 1
            data.group[i] = j
        else
            global j = j + 1
            data.group[i] = j
        end
    end
end

# Group Means and Standard Error
mean_sif740, mean_sif740_sim, stderr_sif740  = Float64[], Float64[], Float64[]
mean_sif757, mean_sif757_sim, stderr_sif757  = Float64[], Float64[], Float64[]
mean_sif771, mean_sif771_sim, stderr_sif771  = Float64[], Float64[], Float64[]
mean_ref757, mean_ref757_sim, stderr_ref757  = Float64[], Float64[], Float64[]
mean_ref771, mean_ref771_sim, stderr_ref771  = Float64[], Float64[], Float64[]
mean_raa, mean_vza  = Float64[], Float64[]
for i in 1:length(unique(data.group))
    if length(data[data[:group] .== [i], :].group) > 9
        m_sif740 = mean(data[data[:group] .== [i], :].sif740)
        m_sif740_sim = mean(data[data[:group] .== [i], :].sif740_sim)
        err_sif740 = std(data[data[:group] .== [i], :].sif740) / sqrt(nrow(data[data[:group] .== [i], :]))
        append!(mean_sif740, m_sif740)
        append!(mean_sif740_sim, m_sif740_sim)
        append!(stderr_sif740, err_sif740)

        m_sif757 = mean(data[data[:group] .== [i], :].sif757)
        m_sif757_sim = mean(data[data[:group] .== [i], :].sif757_sim)
        err_sif757 = std(data[data[:group] .== [i], :].sif757) / sqrt(nrow(data[data[:group] .== [i], :]))
        append!(mean_sif757, m_sif757)
        append!(mean_sif757_sim, m_sif757_sim)
        append!(stderr_sif757, err_sif757)

        m_sif771 = mean(data[data[:group] .== [i], :].sif771)
        m_sif771_sim = mean(data[data[:group] .== [i], :].sif771_sim)
        err_sif771 = std(data[data[:group] .== [i], :].sif771) / sqrt(nrow(data[data[:group] .== [i], :]))
        append!(mean_sif771, m_sif771)
        append!(mean_sif771_sim, m_sif771_sim)
        append!(stderr_sif771, err_sif771)

        m_ref757 = mean(data[data[:group] .== [i], :].rad757)
        m_ref757_sim = mean(data[data[:group] .== [i], :].ref757_sim)
        err_ref757 = std(data[data[:group] .== [i], :].rad757) / sqrt(nrow(data[data[:group] .== [i], :]))
        append!(mean_ref757, m_ref757)
        append!(mean_ref757_sim, m_ref757_sim)
        append!(stderr_ref757, err_ref757)

        m_ref771 = mean(data[data[:group] .== [i], :].rad771)
        m_ref771_sim = mean(data[data[:group] .== [i], :].ref771_sim)
        err_ref771 = std(data[data[:group] .== [i], :].rad771) / sqrt(nrow(data[data[:group] .== [i], :]))
        append!(mean_ref771, m_ref771)
        append!(mean_ref771_sim, m_ref771_sim)
        append!(stderr_ref771, err_ref771)

        m_raa = mean(data[data[:group] .== [i], :].raa)
        m_vza = mean(data[data[:group] .== [i], :].vza)
        append!(mean_raa, m_raa)
        append!(mean_vza, m_vza)

    else
        println("Group ", i, " excluded with an N of ", length(data[data[:group] .== [i], :].group))
    end
end
println("Number of groups included: ", length(mean_sif740))

# Mean aboslute error
mae_sif740 = round(mae(mean_sif740, mean_sif740_sim), digits = 2)
mae_sif757 = round(mae(mean_sif757, mean_sif757_sim), digits = 2)
mae_sif771 = round(mae(mean_sif771, mean_sif771_sim), digits = 2)
mae_ref757 = round(mae(mean_ref757, mean_ref757_sim), digits = 2)
mae_ref771 = round(mae(mean_ref771, mean_ref771_sim), digits = 2)

# Mean Scatter SIF 740, 757, 771
scatter(mean_sif740_sim, mean_sif740, xlabel = "Mean CliMA SIF (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 SIF (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
        title = "ecostress us syv \nLAI = Copernicus, Cab = 40, Radiation = ERA5", titlefontsize = 13, legend = :topleft, label = "740 nm", framestyle = :box, yerror = stderr_sif740)
scatter!(mean_sif757_sim, mean_sif757, label = "757 nm", reg = true, linewidth = 2, yerror = stderr_sif757)
scatter!(mean_sif771_sim, mean_sif771, label = "771 nm", reg = true, linewidth = 2, yerror = stderr_sif771)
df740, df757, df771 = DataFrame(X = mean_sif740_sim, Y = mean_sif740), DataFrame(X = mean_sif757_sim, Y = mean_sif757), DataFrame(X = mean_sif771_sim, Y = mean_sif771)
reg740, reg757, reg771 = lm(@formula(Y ~ X), df740), lm(@formula(Y ~ X), df757), lm(@formula(Y ~ X), df771)
slope740, slope757, slope717 = round(coef(reg740)[2], digits = 2), round(coef(reg757)[2], digits = 2), round(coef(reg771)[2], digits = 2)
intercept740, intercept757, intercept771 = round(coef(reg740)[1]), round(coef(reg757)[1]), round(coef(reg771)[1])
pval740, pval757, pval771 = round(coeftable(reg740).cols[4][2], digits = 3), round(coeftable(reg757).cols[4][2], digits = 3), round(coeftable(reg771).cols[4][2], digits = 3)
if pval757 == 0
    pval757 = "p-value ≤ 0.001"
    pval771 = "p-value ≤ 0.001"
    pval740 = "p-value ≤ 0.001"
else
    pval757 = "p-value = $pval740"
    pval771 = "p-value = $pval757"
    pval740 = "p-value = $pval771"
end
annotate!(3.40, 0.7, text("R² = $(round(r2(reg740), digits = 2)) \nMAE = $mae_sif740 \n$pval740", :left, 10))
annotate!(2.75, 0.7, text("R² = $(round(r2(reg757), digits = 2)) \nMAE = $mae_sif757 \n$pval757", :left, 10))
annotate!(1.75, 0.7, text("R² = $(round(r2(reg771), digits = 2)) \nMAE = $mae_sif771 \n$pval771", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF_OCO3_CliMA_Scatter_LAI-COP_Rad_Mean_ecostress_us_syv.pdf")

# Mean Scatter REF757 and 771
scatter(mean_ref757_sim, mean_ref757, xlabel = "Mean CliMA Reflected Radiance (W m⁻² µm⁻¹ sr⁻¹)", ylabel = "Mean OCO3 Reflected Radiance (W m⁻² µm⁻¹ sr⁻¹)", reg = true, linewidth = 2,
        title = "ecostress us syv \nLAI = Copernicus, Cab = 40, Radiation = ERA5", titlefontsize = 13, legend = :topleft, label = "757 nm", framestyle = :box, yerror = stderr_ref757)
scatter!(mean_ref771_sim, mean_ref771, label = "771 nm", reg = true, linewidth = 2, yerror = stderr_ref771)
df757, df771 = DataFrame(X = mean_ref757_sim, Y = mean_ref757), DataFrame(X = mean_ref771_sim, Y = mean_ref771)
reg757, reg771 = lm(@formula(Y ~ X), df757), lm(@formula(Y ~ X), df771)
slope757, slope717 = round(coef(reg757)[2], digits = 2), round(coef(reg771)[2], digits = 2)
intercept757, intercept771 = round(coef(reg757)[1]), round(coef(reg771)[1])
pval757, pval771 = round(coeftable(reg757).cols[4][2], digits = 3), round(coeftable(reg771).cols[4][2], digits = 3)
if pval757 == 0
    pval757 = "p-value ≤ 0.001"
    pval771 = "p-value ≤ 0.001"
else
    pval757 = "p-value = $pval757"
    pval771 = "p-value = $pval771"
end
annotate!(167.5, 90, text("R² = $(round(r2(reg757), digits = 2)) \nMAE = $mae_ref757 \n$pval757", :left, 10))
annotate!(157.4, 105, text("R² = $(round(r2(reg771), digits = 2)) \nMAE = $mae_ref771 \n$pval771", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/REF_OCO3_CliMA_Scatter_LAI-COP_Rad_Mean_ecostress_us_syv.pdf")

# Polar plot mean OCO3 740 SIF
pyplot()
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(data[:raa]), data[:vza], c = data[:sif740], cmap = :viridis, zorder = 10)
PyPlot.title("OCO-3 SIF 740 nm (SZA 29 - 31.5)")
PyPlot.colorbar(label = "SIF 740 nm (W m⁻² µm⁻¹ sr⁻¹)")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/OCO-3_SIF740_Polar_ecostress_us_syv.pdf")

# Polar plot CliMA 740 SIF
pyplot()
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(data.raa), data.vza, c = data.sif740_sim, cmap = :viridis, zorder = 10)
PyPlot.title("Mean CliMA SIF 740 nm (SZA 29 - 31.5)")
PyPlot.colorbar(label = "SIF 740 nm (W m⁻² µm⁻¹ sr⁻¹)")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/CliMA_SIF740_Polar.pdf")


# Polar plot mean CliMA 757 reflectance
pyplot()
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(mean_raa), mean_vza, c = mean_ref757_sim, cmap = :viridis, zorder = 10)
PyPlot.title("Mean CliMA Reflected Radiance 757 nm (SZA 28.4 - 30.3)")
PyPlot.colorbar(label = "Reflectance (W m⁻² µm⁻¹ sr⁻¹)")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/CliMA_Reflectance_757_Polar.pdf")

# Polar plot mean OCO3 757 reflectance
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(mean_raa), mean_vza, c = mean_ref757, cmap = :viridis, zorder = 10)
PyPlot.title("Mean OCO-3 Reflected Radiance 757 nm (SZA 28.4 - 30.3)")
PyPlot.colorbar(label = "Reflectance (W m⁻² µm⁻¹ sr⁻¹)")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/OCO-3_Reflectance_757_Polar.pdf")

# Polar plot difference 757 reflectance
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(mean_raa), mean_vza, c = (mean_ref757 - mean_ref757_sim), cmap = :viridis, zorder = 10)
PyPlot.title("Difference in Mean OCO-3 and CliMA Reflected Radiance 757 nm (SZA 28.4 - 30.3)")
PyPlot.colorbar(label = "Reflectance (W m⁻² µm⁻¹ sr⁻¹)")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/OCO-3_Reflectance_757_Polar.pdf")

# Polar plot mean CliMA 740 SIF
pyplot()
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(mean_raa), mean_vza, c = mean_sif740_sim, cmap = :viridis, zorder = 10)
PyPlot.title("ecostress us syv \nMean CliMA SIF 740 nm (SZA 29 - 31.5)")
PyPlot.colorbar(label = "SIF 740 nm (W m⁻² µm⁻¹ sr⁻¹)")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/CliMA_SIF740_Polar.pdf")

# Polar plot mean OCO3 740 SIF
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(mean_raa), mean_vza, c = mean_sif740, cmap = :viridis, zorder = 10)
PyPlot.title("ecostress us syv \nMean OCO-3 SIF 740 nm (SZA 29 - 31.5)")
PyPlot.colorbar(label = "SIF 740 nm (W m⁻² µm⁻¹ sr⁻¹)")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/OCO-3_SIF740_Polar.pdf")

# Polar plot difference 740 SIF
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(mean_raa), mean_vza, c = (mean_sif740 - mean_sif740_sim), cmap = :viridis, zorder = 10)
PyPlot.title("Difference in Mean OCO-3 and CliMA SIF 740 nm (SZA 28.4 - 30.3)")
PyPlot.colorbar(label = "SIF 740 nm (W m⁻² µm⁻¹ sr⁻¹)")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/OCO-3_SIF740_Polar.pdf")