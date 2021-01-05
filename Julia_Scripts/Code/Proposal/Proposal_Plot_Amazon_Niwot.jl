code_dir = "C:/Russell/Projects/Geometry/Julia_Scripts/Code/"
include(code_dir*"Derive/derive_spectrum.jl")
include(code_dir*"Derive/group_data_and_stats.jl")

using CSV

amazon_data        = "C:/Russell/Projects/Geometry/R_Scripts/CSV/sif_ATTO_Tower_Manaus_Brazil_(incorrect)_2020-06-26_6500.csv"
amazon_cab         = 60
amazon_week        = NaN
amazon_clumping    = 1
amazon_sif_yield   = 0.5

# Run model
amazon_mat_REF, amazon_mat_SIF, wlf, wl, amazon_data = derive_spectrum(amazon_data, amazon_cab, amazon_week, amazon_clumping, amazon_sif_yield)

# Get mean and mae
amazon_mean, amazon_mae = means_and_errors(dfs)

niwot_mean = CSV.read()
amazon_mae = means_and_errors(dfs)