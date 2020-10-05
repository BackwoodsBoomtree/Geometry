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
    _file = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_ATTO_Tower_Manaus_Brazil_(incorrect)_2020-06-26.csv"
    _data = DataFrame!(CSV.File(_file));

    # initialize canopy radiation module
    angles, can, can_opt, can_rad, in_rad,
    leaves, rt_con, rt_dim, soil, wls = initialize_rt_module(FT);

    # create a copy of incoming radiation
    in_rad_bak = deepcopy(in_rad);
    e_all_dire = sum(in_rad_bak.E_direct  .* wls.dWL) / 1000;
    e_all_diff = sum(in_rad_bak.E_diffuse .* wls.dWL) / 1000;

    # create a matrix to store the spectrum
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

        # change incoming radiation, comment these to see the improvements
        in_rad = deepcopy(in_rad_bak);
        in_rad.E_direct  .*= _data.incoming_direct_era5[i]  / e_all_dire;
        in_rad.E_diffuse .*= _data.incoming_diffuse_era5[i] / e_all_diff;

        for i in 1:20
            leaves[i].Cab = 60;
            fluspect!(leaves[i], wls);
        end

        # re-run the simulations
        canopy_geometry!(can, angles, can_opt, rt_con);
        canopy_matrices!(leaves, can_opt);
        short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
        canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
        SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

        mat_REF[i,:] .= can_rad.alb_obs;
        mat_SIF[i,:] .= can_rad.SIF_obs;
    end

    # 742 nm and 737 nm to 740 nm
    # 757 nm to 757 nm
    # 767 nm and 774.5 nm to 771 nm
    _data.sif740_sim = mat_SIF[:,19] .* 0.4 .+ mat_SIF[:,20] .* 0.6;
    _data.sif757_sim = mat_SIF[:,23];
    _data.sif771_sim = mat_SIF[:,25] .* 0.4667 .+ mat_SIF[:,26] .* 0.5333;

    return mat_REF, mat_SIF, wls.WLF, wls.WL, _data
    #return nothing
end

mat_REF, mat_SIF, wlf, wl, data = derive_spectrum()

# _fig1,_axes = create_canvas(1; figsize=(7.5,3.5));
# _ax1 = _axes[1];
# _cm1 = _ax1.pcolor(mat_REF);
# _ax1.set_title("Reflectance Spectrum");
# _ax1.set_xlabel("Wave Length Index");
# _ax1.set_ylabel("Site Index");
# _fig1.colorbar(_cm1, ax=_ax1, label="Albedo");
# _fig1.set_tight_layout(true);

# _fig2,_axes = create_canvas(2; figsize=(5.5,3.5));
# _ax2 = _axes[1];
# _cm2 = _ax2.pcolor(wlf, eachindex(data.vza), mat_SIF);
# _ax2.set_title("SIF Spectrum");
# _ax2.set_xlabel("Wave Length (nm)");
# _ax2.set_ylabel("Site Index");
# _fig2.colorbar(_cm2, ax=_ax2, label="SIF");
# _fig2.set_tight_layout(true);

# Remove low LAI
data = data[(data[:lai_cop].>4), :]

# Groups
insertcols!(data, ncol(data)+1, :group => NaN)
for i in 1:nrow(data)
    if data[:raa][i] < 45
        data.group[i] = 1
    elseif (data[:raa][i] < 65) && (data[:raa][i] > 45)
        data.group[i] = 2
    elseif (data[:raa][i] < 80) && (data[:raa][i] > 65)
        data.group[i] = 3
    elseif (data[:raa][i] < 135) && (data[:raa][i] >90)
        data.group[i] = 4
    elseif (data[:raa][i] < 157) && (data[:raa][i] > 150) && (data[:vza][i] > 25)
        data.group[i] = 5
    elseif (data[:raa][i] < 171) && (data[:raa][i] > 157) && (data[:vza][i] > 35)
        data.group[i] = 6
    elseif (data[:raa][i] > 174) && (data[:vza][i] > 40)
        data.group[i] = 7
    elseif (data[:raa][i] > 174) && (data[:vza][i] < 40)
        data.group[i] = 8
    elseif (data[:raa][i] < 171) && (data[:raa][i] > 157) && (data[:vza][i] < 35)
        data.group[i] = 9
    end
end

# Group Means and Standard Error
mean_sif740 = Float64[]
mean_sif740_sim = Float64[]
stderr_sif740 = Float64[]
for i in 1:length(unique(data.group))
    if length(data[data[:group] .== [i], :].group) > 9
        m_sif = mean(data[data[:group] .== [i], :].sif740)
        m_sif_sim = mean(data[data[:group] .== [i], :].sif740_sim)
        err_sif = std(data[data[:group] .== [i], :].sif740) / sqrt(nrow(data[data[:group] .== [i], :]))
        append!(mean_sif740, m_sif)
        append!(mean_sif740_sim, m_sif_sim)
        append!(stderr_sif740, err_sif)
    else
        println("Group ", i, " excluded with an N of ", length(data[data[:group] .== [i], :].group))
    end
end
println("Number of groups included: ", length(mean_sif740))

# Mean aboslute error
sif_mae = round(mae(mean_sif740, mean_sif740_sim), digits = 2)

# Mean Scatter
scatter(mean_sif740_sim, mean_sif740, xlabel = "Mean CliMA SIF740", ylabel = "Mean OCO3 SIF740",
        title = "LAI = Copernicus, Cab = 60, Radiation = Default", titlefontsize = 13, legend = false, framestyle = :box, yerror = stderr_sif740)
df = DataFrame(X = mean_sif740_sim, Y = mean_sif740)
reg = lm(@formula(Y ~ X), df)
slope = round(coef(reg)[2], digits = 2)
intercept = round(coef(reg)[1])
pval = round(coeftable(reg).cols[4][2], digits = 3)
if pval == 0
    pval = "p-value ≤ 0.001"
else
    pval = "p-value = $pval"
end
plot!(mean_sif740_sim, predict(reg), w = 3)
annotate!(:topleft, text("R² = $(round(r2(reg), digits = 2)) \nMAE = $sif_mae \n$pval", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF740_OCO3_CliMA_Scatter_LAI-COP_Mean.pdf")

# LAI Scatter Plots
scatter(sif740, data.sif740, xlabel = "CliMA SIF740", ylabel = "OCO3 SIF740", xlims = (2.25, 3), ylims = (0, 3),
        title = "LAI = Copernicus, Cab = 60, Radiation = Default", titlefontsize = 13, legend = false, colorbar = true, framestyle = :box, zcolor = data.lai_cop, m = (:viridis, 0.8),
        colorbar_title = "LAI")
df = DataFrame(X = sif740, Y = data.sif740)
reg = lm(@formula(Y ~ X), df)
slope = round(coef(reg)[2], digits = 2)
intercept = round(coef(reg)[1])
pval = round(coeftable(reg).cols[4][2], digits = 3)
if pval == 0
    pval = "p-value ≤ 0.001"
else
    pval = "p-value = $pval"
end
plot!(sif740, predict(reg), w = 3)
annotate!(2.30, 2.75, text("R² = $(round(r2(reg), digits = 2)) \n$pval", :left, 10))
#annotate!(2.51, 2.75, text("R² = $(round(r2(reg), digits = 2)) \n$pval \ny = $slope * x + $intercept", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF740_OCO3_CliMA_Scatter_LAI-COP.pdf")




_fig3,_axes = create_canvas(3; figsize=(3.5,3.5));
_ax3 = _axes[1];
_ax3.plot(_sif740, data.sif740, "k+");
plot_line_regress(_ax3, _sif740, convert(Array, data.sif740); interval=true);
_ax3.set_xlabel("SIF 740");
_ax3.set_ylabel("SIF 740");
_fig3.set_tight_layout(true);

_fig4,_axes = create_canvas(4; figsize=(3.5,3.5));
_ax4 = _axes[1];
_ax4.plot(_sif757, data.sif757, "k+");
plot_line_regress(_ax4, _sif757, convert(Array, data.sif757); interval=true);
_ax4.set_xlabel("SIF 757");
_ax4.set_ylabel("SIF 757");
_fig4.set_tight_layout(true);

_fig5,_axes = create_canvas(5; figsize=(3.5,3.5));
_ax5 = _axes[1];
_ax5.plot(_sif771, data.sif771, "k+");
plot_line_regress(_ax5, _sif771, convert(Array, data.sif771); interval=true);
_ax5.set_xlabel("SIF 771");
_ax5.set_ylabel("SIF 771");
_fig5.set_tight_layout(true);