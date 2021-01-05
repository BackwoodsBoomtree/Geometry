using Plots
using Plots.PlotMeasures
using Parameters
using Land
using Land.CanopyRT
#@unpack FT, leafbio, canopy, angles, canOpt, canRad,sunRad,soil, wl, wle, wlf = CanopyRT

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
for i in 1:canopy_rt.nlayers
    fluspect!(arrayOfLeaves[i],  wl_set)
end

# Basic Steps for Canopy RT ###
# Set Soil albedo to 0.2
soil = SoilOpti{FT}(wl_set.wl, FT(0.2)*ones(FT, length(wl_set.wl)), FT[0.1], FT(290.0))
angles = SolarAngles{FT}()

# Compute Canopy optical properties dependend on sun-sensor and leaf angle distributions:
compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
# Compute RT matrices with leaf reflectance and transmissions folded in:
compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
# Perform SW radiation transfer:
simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
# Compute outgoing SIF flux (using constant fluorescence efficiency at the chloroplast level)
derive_canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set);

# Test a VZA dependence in the principal plane ####
SIF_FR = Float32[]
SIF_R = Float32[]
reflVIS = Float32[]
reflNIR = Float32[]

# Just running the code over all geometries:
# Set sun SZA to 30 degrees
angles.tts = 30
# Set 0 azimuth (principal plane)
angles.psi = 0
# LAI of 3:
canopy_rt.LAI = 3
# Define VZA
VZA = collect(-89.5:0.5:89.5)

for VZA_ in VZA
    angles.tto = VZA_
    compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
    compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
    simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
    computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);
    # Handpicked indices in
    push!(reflVIS, canRad_rt.alb_obs[ind_red])
    push!(reflNIR, canRad_rt.alb_obs[ind_NIR])
    push!(SIF_R , canRad_rt.SIF_obs[ind_wlf_R])
    push!(SIF_FR, canRad_rt.SIF_obs[ind_wlf_FR ])
end

# BRDF sampling
# By going through viewing and azimuth angles, we can construct a full BRDF for reflectance and SIF emissions at different wavelengths
reflVIS = Float32[]
reflNIR = Float32[]
SIF_FR = Float32[]
SIF_R  = Float32[]
angles.tts = 30
angles.psi = 0
canopy_rt.LAI = 3.22
for psi = 0:360
    angles.psi = psi
    for VZA = 0:1:85
        angles.tto = VZA

        compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
        compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
        simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
        computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);
        push!(reflVIS, canRad_rt.alb_obs[28])
        push!(reflNIR, canRad_rt.alb_obs[52])
        push!(SIF_R , canRad_rt.SIF_obs[8])
        push!(SIF_FR, canRad_rt.SIF_obs[20])
    end
end

A = reshape(reflNIR, ( 86,361));
B = reshape(reflVIS, ( 86,361));
SIFFER = reshape(SIF_R, ( 86,361));
SIFFER_FR = reshape(SIF_FR, ( 86,361));

# 3D Plot
x, y = 0:180, 0:85
#z = A
z = A[1:86,1:181]
z2 = SIFFER_FR[1:86,1:181]/6

pyplot()
surface(ylabel = "Viewing Zenith Angle (°)", xlabel = "Relative Azimuth Angle (°)", zlabel = "Emission or Reflectance",
        zguidefontrotation = 90, yguidefontrotation = -46, xguidefontrotation = 18,
        title = "Solar Zenith Angle = 30°", size = (800,600), camera=(-30,30))
p1 = surface!(x,y,z, colorbar = false)
p2 = surface!(x,y,z2, c = :viridis, colorbar_title = "SIF/6", colorbar = true)

x = collect(0:360)
y = collect(0:85)
z = fill(30, (86, 361))

pyplot()
surface(ylabel = "Viewing Zenith Angle (°)", xlabel = "Relative Azimuth Angle (°)", zlabel = "Emission or Reflectance",
        title = "Solar Zenith Angle = 30°", size = (800,600))
surface!(x,y,z, fill_z = sif_757)
surface!(fill_z = sif_757, c = :viridis, colorbar_title = "asdf", colorbar = true)