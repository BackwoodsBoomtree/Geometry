include("derive_spectrum_sensitivity.jl")
include("Plot_Sensitivity_Curves.jl")

model_variable = "rad_ratio"
range = [1,7.0]
by = 2.0
filename = "C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Sensitivity/Radiation_Ratio.pdf"

mat_REF, mat_SIF, wlf, wl, df = derive_spectrum_sensitivity(model_variable, range, by)

p_sensitivty_curves = sensitivity_curves(model_variable, df, mat_SIF, mat_REF, wlf, wl)

savefig(p_sensitivty_curves, filename)