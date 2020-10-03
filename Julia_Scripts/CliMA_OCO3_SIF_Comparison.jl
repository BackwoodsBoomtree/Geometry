using CSV
using Plots
using Plots.PlotMeasures
using Parameters
using Land
using Land.CanopyRT
# using CanopyLayers
using GLM
using Statistics
using StatsPlots
using DataFrames

oco3_data = CSV.read("C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_ATTO_Tower_Manaus_Brazil_(incorrect)_2020-06-26.csv")

const FT = Float32

wl_set = create_wl_para_set(FT)
leaf = create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF);
canopy_rt = Canopy4RT{FT, 20, 3.0}()
canRad_rt = CanopyRadiation{FT, wl_set.nwl, wl_set.nWlF, length(canopy_rt.litab), length(canopy_rt.lazitab), canopy_rt.nlayers}()
canOpt_rt = create_canopy_optical(FT, wl_set.nwl, canopy_rt.nlayers, length(canopy_rt.lazitab), length(canopy_rt.litab); using_marray=false)
sunRad_rt = create_incoming_radiation(FT, wl_set.swl);

# Run Fluspect:
fluspect!(leaf, wl_set);

# Define a few wavelengths:
wl_blue = 450.0;
wl_red = 600.0;
wl_FarRed = 740.0;
wl_Red = 685.0;
ind_wle_blue  = argmin(abs.(wl_set.wle .-wl_blue));
ind_wle_red = argmin(abs.(wl_set.wle .-wl_red));
ind_wlf_FR  = argmin(abs.(wl_set.wlf .-wl_FarRed));
ind_wlf_R  = argmin(abs.(wl_set.wlf .-wl_Red));
ind_red = argmin(abs.(wl_set.wl .-wl_Red));
ind_NIR = argmin(abs.(wl_set.wl .-800));

arrayOfLeaves = [create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF) for i in 1:canopy_rt.nlayers]
for cab = 1:length(arrayOfLeaves)
    arrayOfLeaves[cab].Cab = 60
end

for i in 1:canopy_rt.nlayers
    fluspect!(arrayOfLeaves[i],  wl_set)
end

# Basic Steps for Canopy RT ###
# Set Soil albedo to 0.2
soil = SoilOpti{FT}(wl_set.wl, FT(0.2)*ones(FT, length(wl_set.wl)), FT[0.1], FT(290.0))
angles = SolarAngles{FT}()

#### Sun-sensor ####
reflVIS = Float32[]
reflNIR = Float32[]
SIF_FR = Float32[]
SIF_R  = Float32[]

# Run code over OCO geometries
SZA = oco3_data[:sza]
RAA = oco3_data[:raa]
VZA = oco3_data[:vza]

# LAI
LAI = oco3_data[:lai]

for i = 1:length(VZA)
    angles.tts = SZA[i]
    angles.psi = RAA[i]
    angles.tto = VZA[i]
    canopy_rt.LAI = LAI[i]
    # canopy_rt.LAI = 5

    # for j = 1:length(arrayOfLeaves)
    #     arrayOfLeaves[j].Cab = arrayOfLeaves[j].Cab * EVInorm[i]
    # end

    compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
    compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
    simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
    computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);
    push!(reflVIS, canRad_rt.alb_obs[28])
    push!(reflNIR, canRad_rt.alb_obs[52])
    push!(SIF_R  , (canRad_rt.SIF_obs[19] + canRad_rt.SIF_obs[20]) / 2)
    push!(SIF_FR, canRad_rt.SIF_obs[26])
end

# Plot Red
plot(framestyle = :box, xlabel = "Viewing Zenith Angle (degrees)", ylabel = "Reflectance or Emission", legend = :topleft)
plot!(VZA, reflVIS, label = "Red Reflectance", line = 2)
plot!(VZA, SIF_R/30, label = "Red SIF (/30)", line = 2)

# Plot Near Infrared:
plot(framestyle = :box, xlabel = "Viewing Zenith Angle (degrees)", ylabel = "Reflectance or Emission", legend = :topleft)
plot!(VZA, reflNIR, label = "NIR Reflectance", line = 2)
plot!(VZA, SIF_FR/6, label = "Far Red SIF (/6)", line = 2)

# Need to convert reclectance/emission vectors to 2D arrays
reflNIR_matrix = Array{Float64}(undef, length(VZA), length(RAA))
replace!(reflNIR_matrix, 0 => NaN)
SIF_FR_matrix = Array{Float64}(undef, length(VZA), length(RAA))
replace!(SIF_FR_matrix, 0 => NaN)
for i = 1:length(reflNIR)
    reflNIR_matrix[i,i] = reflNIR[i]
    SIF_FR_matrix[i,i] = SIF_FR[i]
end

# Polar plot NIR CliMA
pyplot()
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(RAA), VZA, c = reflNIR, cmap = :viridis, zorder = 10)
PyPlot.title("CliMA - NIR Reflectance (SZA 28.4 - 30.3)")
PyPlot.colorbar(label = "NIR Reflectance")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Radiance_NIR_CliMA_Polar.pdf")

# Polar plot SIF740 CliMA
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(RAA), VZA, c = SIF_FR, cmap = :viridis, zorder = 10)
PyPlot.title("CliMA - SIF740 (SZA 28.4 - 30.3)")
PyPlot.colorbar(label = "SIF740")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF740_CliMA_Polar.pdf")

# Polar plot Radiance OCO
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(RAA), VZA, c = oco3_data[:sif771_R], cmap = :viridis, zorder = 10)
PyPlot.title("OCO3 - 771nm Continnum Level Radiance (SZA 28.4 - 30.3)")
PyPlot.colorbar(label = "771nm Continnum Level Radiance")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Radiance_NIR_OCO3_Polar.pdf")

# Polar plot SIF740 OCO
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(true)
hm = PyPlot.scatter(deg2rad.(RAA), VZA, c = oco3_data[:sif740], cmap = :viridis, zorder = 10)
PyPlot.title("OCO3 - SIF740 (SZA 28.4 - 30.3)")
PyPlot.colorbar(label = "SIF740")
PyPlot.gcf()
PyPlot.savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF740_OCO3_Polar.pdf")

# Scatter Plots
scatter(SIF_FR, oco3_data[:sif740], xlabel = "CliMA SIF740", ylabel = "OCO3 SIF740", xlims = (2.5, 3.1), ylims = (0, 3.1),
        title = "LAI = 1 km Copernicus, Cab = 60", titlefontsize = 13, legend = false, colorbar = true, framestyle = :box, zcolor = LAI, m = (:viridis, 0.8),
        colorbar_title = "LAI")
data = DataFrame(X = SIF_FR, Y = oco3_data[:sif740])
reg = lm(@formula(Y ~ X), data)
slope = round(coef(reg)[2], digits = 2)
intercept = round(coef(reg)[1])
pval = round(coeftable(reg).cols[4][2], digits = 3)
if pval == 0
    pval = "p-value ≤ 0.001"
else
    pval = "p-value = $pval"
end
plot!(SIF_FR, predict(reg), w = 3)
annotate!(2.51, 2.75, text("R² = $(round(r2(reg), digits = 2)) \n$pval", :left, 10))
#annotate!(2.51, 2.75, text("R² = $(round(r2(reg), digits = 2)) \n$pval \ny = $slope * x + $intercept", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/SIF740_OCO3_CliMA_Scatter_LAI.pdf")

scatter(reflNIR, oco3_data[:rad771], xlabel = "CliMA NIR Reflectance", ylabel = "OCO3 771nm Continnum Level Radiance", xlims = (0.40, 0.52), ylims = (0, 150),
        title = "LAI = 1 km Copernicus, Cab = 60", titlefontsize = 13, legend = false, colorbar = true, framestyle = :box, zcolor = LAI, m = (:viridis, 0.8),
        colorbar_title = "LAI")
data = DataFrame(X = reflNIR, Y = oco3_data[:rad771])
reg = lm(@formula(Y ~ X), data)
slope = round(coef(reg)[2], digits = 2)
intercept = round(coef(reg)[1])
pval = round(coeftable(reg).cols[4][2], digits = 3)
if pval == 0
    pval = "p-value ≤ 0.001"
else
    pval = "p-value = $pval"
end
plot!(reflNIR, predict(reg), w = 3)
annotate!(0.405, 138, text("R² = $(round(r2(reg), digits = 2)) \n$pval", :left, 10))
#annotate!(0.405, 138, text("R² = $(round(r2(reg), digits = 2)) \n$pval \ny = $slope * x + $intercept", :left, 10))
savefig("C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Radiance_NIR_OCO3_CliMA_Scatter_LAI.pdf")

# Min and Max SZA from oco3 data
szamin = minimum(oco3_data[:sza])
szamax = maximum(oco3_data[:sza])
string("SZA min/max: ", szamin, " ", szamax)

# Min and Max VZA from oco3 data
vzamin = minimum(oco3_data[:vza])
vzamax = maximum(oco3_data[:vza])
string("VZA min/max: ", vzamin, " ", vzamax)

# Min and Max PA from oco3 data
pamin = minimum(oco3_data[:pa])
pamax = maximum(oco3_data[:pa])
string("PA min/max: ", pamin, " ", pamax)

# Min and Max RAA from oco3 data
raamin = minimum(oco3_data[:raa])
raamax = maximum(oco3_data[:raa])
string("RAA min/max: ", raamin, " ", raamax)

# Min and Max VAA from oco3 data
vaamin = minimum(oco3_data[:vaa])
vaamax = maximum(oco3_data[:vaa])
string("VAA min/max: ", vaamin, " ", vaamax)

