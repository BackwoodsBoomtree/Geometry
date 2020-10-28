using DataFrames
using NCDatasets

function compute_relative_angle(df)
    # Calculate Relative Azimuth Angle
    df.RAz = broadcast(abs, (df.SAz - df.VAz))
    for a in 1:length(df.RAz)
        if df.RAz[a] > 180
            df.RAz[a] = abs(df.RAz[a] - 360)
        end
    end
end

function compute_phase_angle(df)
    ## Phase angles are set to negative values if the observational azimuth angle is bigger than the solar azimuth angle
    ## (negative phase angle means the sun is to the right of the satellite)
    compute_relative_angle(df)
    phase = Int64[]
    df.PA = 0.0
    for i in 1:length(df.SAz)
        if df.VAz[i] > df.SAz[i]
            df.PA[i] = -1
        elseif df.VAz[i] < df.SAz[i]
            df.PA[i] = 1
        end
    end
    for i in 1:length(df.PA)
        df.PA[i] = (acos(cosd(df.SZA[i]) * cosd(df.VZA[i]) + sind(df.VZA[i]) * sind(df.SZA[i]) * cosd(df.RAz[i])) * 180 / π) * df.PA[i]
    end
    return(df)
end

function build_oco3_data(file_list)
    # Define which variables to include in output
    var_groups = ["Cloud", "Geolocation", "Metadata", "Meteo", "Offset", "Science"] # Groups from NC file
    var_cloud = ["cloud_flag_abp"]
    var_geoloc = ["longitude", "latitude"]
    var_meta = ["MeasurementMode", "OrbitId"]
    var_science = ["IGBP_index", "sounding_land_fraction", "continuum_radiance_757nm", "continuum_radiance_771nm", "SIF_757nm", "SIF_771nm", "SIF_Relative_757nm", "SIF_Relative_771nm"]
    var_main = ["Quality_Flag", "SAz", "SZA", "VAz", "VZA", "SIF_740nm", "Delta_Time"] # Main variables; not grouped

    df = DataFrame()
    for f in 1:length(list_files)
        df_temp = DataFrame()
        ncd = NCDataset(list_files[f])
        if f == 1
            for v in 1:length(var_main)
                df.var = ncd[var_main[v]][:,:]
                rename!(df, Dict(:var => Symbol.(var_main[v])))
            end
            for c in 1:length(var_cloud)
                df.var = ncd.group["Cloud"][var_cloud[c]][:,:]
                rename!(df, Dict(:var => Symbol.(var_cloud[c])))
            end
            for g in 1:length(var_geoloc)
                df.var = ncd.group["Geolocation"][var_geoloc[g]][:,:]
                rename!(df, Dict(:var => Symbol.(var_geoloc[g])))
            end
            for m in 1:length(var_meta)
                df.var = ncd.group["Metadata"][var_meta[m]][:,:]
                rename!(df, Dict(:var => Symbol.(var_meta[m])))
            end
            for s in 1:length(var_science)
                df.var = ncd.group["Science"][var_science[s]][:,:]
                rename!(df, Dict(:var => Symbol.(var_science[s])))
            end
        else
            for v in 1:length(var_main)
                df_temp.var = ncd[var_main[v]][:,:]
                rename!(df_temp, Dict(:var => Symbol.(var_main[v])))
            end
            for c in 1:length(var_cloud)
                df_temp.var = ncd.group["Cloud"][var_cloud[c]][:,:]
                rename!(df_temp, Dict(:var => Symbol.(var_cloud[c])))
            end
            for g in 1:length(var_geoloc)
                df_temp.var = ncd.group["Geolocation"][var_geoloc[g]][:,:]
                rename!(df_temp, Dict(:var => Symbol.(var_geoloc[g])))
            end
            for m in 1:length(var_meta)
                df_temp.var = ncd.group["Metadata"][var_meta[m]][:,:]
                rename!(df_temp, Dict(:var => Symbol.(var_meta[m])))
            end
            for s in 1:length(var_science)
                df_temp.var = ncd.group["Science"][var_science[s]][:,:]
                rename!(df_temp, Dict(:var => Symbol.(var_science[s])))
            end
        end
        if f != 1
            append!(df, df_temp)
        end
        print(list_files[f])
    end
    compute_phase_angle(df)
    return(df)
end

function build_era5_data(df, era_location)
    
end