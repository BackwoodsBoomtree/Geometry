code_dir = "C:/Russell/Projects/Geometry/Julia_Scripts/Code/"
include(code_dir*"extract_GOSAT_OCO_data.jl")

file_list_oco2  = ["C:/Russell/Projects/Geometry/Data/oco2/oco2_LtSIF_190104_B10206r_200729173048s.nc4"]
file_list_gosat = ["C:/Russell/Projects/Geometry/Data/gosat/gosat_LtSIF_100301_20200624t080218z.nc4"]

# List of variables from each group in the nc file.
vars_cloud   = ["cloud_flag_abp"]
vars_geoloc  = ["longitude", "latitude"]
vars_meta    = ["MeasurementMode", "OrbitId"]
vars_meteo   = []
vars_science = ["IGBP_index", "sounding_land_fraction", "continuum_radiance_757nm", "continuum_radiance_771nm", "SIF_757nm", "SIF_771nm", "SIF_Relative_757nm", "SIF_Relative_771nm"]
vars_main    = ["Quality_Flag", "SAz", "SZA", "VAz", "VZA", "SIF_740nm", "Delta_Time"] # Variables not grouped

df_oco2  = extract_data(file_list_oco2, vars_cloud, vars_geoloc, vars_meta, vars_meteo, vars_science, vars_main);
df_gosat = extract_data(file_list_gosat, vars_cloud, vars_geoloc, vars_meta, vars_meteo, vars_science, vars_main);