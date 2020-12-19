include("derive_spectrum_sensitivity.jl")
include("Plot_Sensitivity_Curves.jl")

model_variable = "Ω"
range = [0.1,1]
by = 0.1

mat_REF, mat_SIF, wlf, wl, df = derive_spectrum_sensitivity(model_variable, range, by)

p_sensitivty_curves = sensitivity_curves(model_variable, df, mat_SIF, mat_REF, wlf, wl)

