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

# show leaf Chlorophyll content:
@show leaf.Cab

# Direct and diffuse irradiance
plot(title = "Solar Irradiance at Top of Canopy", framestyle = :box, xlabel = "Wavelength (nm)", ylabel = "Irradiance (W/m²)")
plot!(sunRad_rt.wl, sunRad_rt.E_direct, label = "Direct")
plot!(sunRad_rt.wl, sunRad_rt.E_diffuse, label = "Diffuse")

# Run Fluspect:
fluspect!(leaf, wl_set);

# Leaf Reflectance
plot(wl_set.wl, leaf.ρ_SW, xlabel = "Wavelength", ylabel = "Leaf Reflectance", framestyle = :box, leg = false)

# Fluorescence excitation matrices
contourf(wl_set.wle, wl_set.wlf, leaf.Mb, right_margin = 10px,
         xlabel = "Excitation wavelength (nm)",
         ylabel = "Emission wavelength (nm)",
         title = "Fluorescence backward (refl) emission (Mb)")

contourf(wl_set.wle, wl_set.wlf, leaf.Mf, right_margin = 10px,
         xlabel = "Excitation wavelength (nm)",
         ylabel = "Emission wavelength (nm)",
         title = "Fluorescence forward (transmission) emission (Mf)")
         
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

plot(framestyle = :box, size = (800,400), xlabel = "Wavelength (nm)", ylabel = "Fluorescence")
plot!(wl_set.wlf,leaf.Mf[:,ind_wle_blue], color = :black, line = 2 , label = "Forward SIF, excited at $wl_blue nm")
plot!(wl_set.wlf,leaf.Mb[:,ind_wle_blue], color = :orange, line = 2 , label = "Backward SIF, excited at $wl_blue nm")
plot!(wl_set.wlf,leaf.Mf[:,ind_wle_red], color = :black, line = (:dash, 2), label = "Forward SIF, excited at $wl_red nm")
plot!(wl_set.wlf,leaf.Mb[:,ind_wle_red], color = :orange, line = (:dash, 2), label = "Backward SIF, excited at $wl_red nm")

plot(framestyle = :box, size = (800,400), xlabel = "Absorbed Wavelength (nm)", ylabel = "Fluorescence", ylims = (0, 0.005))
plot!(wl_set.wle,leaf.Mf[ind_wlf_FR,:], color = :black, line = 2 , label="Forward SIF, emitted at $wl_FarRed nm")
plot!(wl_set.wle,leaf.Mb[ind_wlf_FR,:],color = :orange, line = 2 , label="Backward SIF, emitted at $wl_FarRed nm")
plot!(wl_set.wle,leaf.Mf[ind_wlf_R,:], color = :black, line = (:dash, 2), label="Forward SIF, emitted at $wl_Red nm")
plot!(wl_set.wle,leaf.Mb[ind_wlf_R,:], color = :orange, line = (:dash, 2), label="Backward SIF, emitted at $wl_Red nm")

# Leaf reflectance and transmission
# Let's create a leaf with a different Cab and Cw (water) content
# Try changing other pigment contents, plot leaf reflectance and transmissions and explain where (spectrally) and why reflectance and transmission changes:
leaf_2 = create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF);
leaf_2.Cab = 80
leaf_2.Cw = 0.012

# Run Fluspect:
fluspect!(leaf_2, wl_set); 

plot(framestyle = :box, size = (800,400), xlabel = "Wavelength (nm)", ylabel = "Transmission or Reflectance", legend = :right)
plot!(wl_set.wl,1 .-leaf.τ_SW, color = :black, line = 2 , label = "Leaf Transmission")
plot!(wl_set.wl,leaf.ρ_SW, color = :orange, line = 2 , label = "Leaf Reflectance" )
plot!(wl_set.wl,1 .-leaf_2.τ_SW, color = :black, line = (:dash, 2), label = "Leaf ##2 Transmission")
plot!(wl_set.wl,leaf_2.ρ_SW, color = :orange, line = (:dash, 2), label = "Leaf ##2 Reflectance" )

# Moving from the leaf to the entire canopy ###
# This is to be changed later but at the moment, we need to generate an Array of leaves, basically for each layer of the canopy
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
angles.tts = 29
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

# Plot Red
plot(framestyle = :box, xlabel = "Viewing Zenith Angle (degrees)", ylabel = "Reflectance or Emission", legend = :topleft)
plot!(VZA, reflVIS, label = "Red Reflectance", line = 2)
plot!(VZA, SIF_R/30, label = "Red SIF (/30)", line = 2)

# Plot Near Infrared:
plot(framestyle = :box, xlabel = "Viewing Zenith Angle (degrees)", ylabel = "Reflectance or Emission", legend = :topleft)
plot!(VZA, reflNIR, label = "NIR Reflectance", line = 2)
plot!(VZA, SIF_FR/6, label = "Far Red SIF (/6)", line = 2)

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

A = reshape(reflNIR, (86,361));
B = reshape(reflVIS, (86,361));
SIFFER = reshape(SIF_R, (86,361));
SIFFER_FR = reshape(SIF_FR, (86,361));

pyplot()
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(false)
hm = PyPlot.contourf(deg2rad.(collect((0:360))),collect(0:1:85),  A,  cmap=:viridis)
PyPlot.title("NIR reflectance BRDF")
PyPlot.colorbar()
PyPlot.gcf()

PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(false)
hm = PyPlot.contourf(deg2rad.(collect((0:360))),collect(0:1:85),  B,  cmap=:viridis)
PyPlot.title("VIS reflectance BRDF")
PyPlot.colorbar()
PyPlot.gcf()

PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(false)
hm = PyPlot.contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER, cmap=:viridis)
PyPlot.title("Red SIF emission BRDF")
PyPlot.colorbar()
PyPlot.gcf()

PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(false)
hm = PyPlot.contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER_FR, cmap=:viridis)
PyPlot.title("Far Red SIF emission BRDF")
PyPlot.colorbar()
PyPlot.gcf()