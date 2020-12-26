include("derive_spectrum_sensitivity.jl")
include("Plot_Sensitivity_Curves.jl")

model_variable = "rad_ratio"
range = [1.0,7.0] # Float (nLayer is Int)
by = 2
filename = "C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Sensitivity/Radiation_Ratio_Soil_1-7.pdf"

mat_rad_obs, mat_ref, mat_alb_obs, mat_SIF_obs, mat_SIF_obs_sunlit, mat_SIF_obs_shaded, 
mat_SIF_obs_scattered, mat_SIF_obs_soil, wlf, wl, df = derive_spectrum_sensitivity(model_variable, range, by)

# Change mat_SIF here for sunlit, shaded, scattered, and soil if interested
p_sensitivty_curves = sensitivity_curves(model_variable, df, mat_SIF_obs_soil, mat_alb_obs, wlf, wl)

savefig(p_sensitivty_curves, filename)

### In addition to the model variables, the following can be run ###
#
# "rad_ratio"     - Compute ratios of direct:diffuse light.
#                   Can use mat_SIF_sunlit, mat_SIF_shaded, mat_SIF_scattered, or mat_SIF_soil
#
# "cab_gradient"  - Run gradients of Cab spread out among the layers.
#                   by argument is not used here as range will be divided among the nLayer
#
# "cab_Yang_2017" - Run gradients of Cab from Yang et al 2017.
#                   range and by arguments are not used