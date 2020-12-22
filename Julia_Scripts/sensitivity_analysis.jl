include("derive_spectrum_sensitivity.jl")
include("Plot_Sensitivity_Curves.jl")

model_variable = "rad_ratio"
range = [1.0,7.0] # Float
by = 2
filename = "C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Sensitivity/Radiation_Ratio_Soil.pdf"

mat_REF, mat_SIF, mat_SIF_sunlit, mat_SIF_shaded, 
mat_SIF_scattered, mat_SIF_soil, wlf, wl, df = derive_spectrum_sensitivity(model_variable, range, by)

# Change mat_SIF here for sunlit, shaded, scattered, and soil if interested
p_sensitivty_curves = sensitivity_curves(model_variable, df, mat_SIF_soil, mat_REF, wlf, wl)

savefig(p_sensitivty_curves, filename)

### In addition to the model variables, the following can be run ###
# "rad_ratio" - Will compute ratios of direct:diffuse light.
#               Can use mat_SIF_sunlit, mat_SIF_shaded, mat_SIF_scattered, or mat_SIF_soil
# "cab_gradient - Will run gradients of Cab spread out among the layers