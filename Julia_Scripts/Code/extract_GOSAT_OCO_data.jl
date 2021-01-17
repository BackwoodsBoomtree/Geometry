code_dir = "C:/Russell/Projects/Geometry/Julia_Scripts/Code/"
include(code_dir*"compute_angles.jl")

using DataFrames
using NCDatasets

# Note: Not meant to be used with data from Offset group

function extract_data(file_list, vars_cloud, vars_geoloc, vars_meta, vars_meteo, vars_science, vars_main)

    df = DataFrame()
    
    for f in 1:length(file_list)
        
        df_temp = DataFrame()
        ncd = NCDataset(file_list[f]);

        for v in 1:length(vars_cloud)
            if haskey(ncd.group["Cloud"], vars_cloud[v])
                if size(ncd.group["Cloud"][vars_cloud[v]][:,:])[1] == 2                        # For GOSAT data (2 polarizations)
                    df[Symbol(vars_cloud[v], "_P")] = ncd.group["Cloud"][vars_cloud[v]][1,:]
                    df[Symbol(vars_cloud[v], "_S")] = ncd.group["Cloud"][vars_cloud[v]][2,:]
                else
                    df[Symbol(vars_cloud[v])] = ncd.group["Cloud"][vars_cloud[v]][:,:]       # For OCO data
                end
            else
                println("The group Cloud does not contain the variable: ", vars_cloud[v])
            end
        end
        for v in 1:length(vars_geoloc)
            if haskey(ncd.group["Geolocation"], vars_geoloc[v])
                df[Symbol(vars_geoloc[v])] = ncd.group["Geolocation"][vars_geoloc[v]][:,:]
            else
                println("The group Geolocation does not contain the variable: ", vars_geoloc[v])
            end
        end
        for v in 1:length(vars_meta)
            if haskey(ncd.group["Metadata"], vars_meta[v])
                df[Symbol(vars_meta[v])] = ncd.group["Metadata"][vars_meta[v]][:,:]
            else
                println("The group Metadata does not contain the variable: ", vars_meta[v])
            end
        end
        for v in 1:length(vars_meteo)
            if haskey(ncd.group["Meteo"], vars_meteo[v])
                df[Symbol(vars_meteo[v])] = ncd.group["Meteo"][vars_meteo[v]][:,:]
            else
                println("The group Meteo does not contain the variable: ", vars_meteo[v])
            end
        end
        for v in 1:length(vars_science)
            if haskey(ncd.group["Science"], vars_science[v])
                if size(ncd.group["Science"][vars_science[v]][:,:])[1] == 2                        # For GOSAT data (2 polarizations)
                    df[Symbol(vars_science[v], "_P")] = ncd.group["Science"][vars_science[v]][1,:]
                    df[Symbol(vars_science[v], "_S")] = ncd.group["Science"][vars_science[v]][2,:]
                else
                    df[Symbol(vars_science[v])] = ncd.group["Science"][vars_science[v]][:,:]       # For OCO data
                end
            else
                println("The group Science does not contain the variable: ", vars_science[v])
            end
        end
        for v in 1:length(vars_main)
            if haskey(ncd, vars_main[v])
                if size(ncd[vars_main[v]][:,:])[1] == 2                        # For GOSAT data (2 polarizations)
                    df[Symbol(vars_main[v], "_P")] = ncd[vars_main[v]][1,:]
                    df[Symbol(vars_main[v], "_S")] = ncd[vars_main[v]][2,:]
                else
                    df[Symbol(vars_main[v])] = ncd[vars_main[v]][:,:]          # For OCO data
                end
            else
                println("The NC file does not contain the ungrouped variable: ", vars_main[v])
            end
        end

        append!(df, df_temp)
        println("Added to the output dataframe: ", file_list[f])

    end
        # Add relative azimuth and phase angles
    compute_relative_angle(df, "SAz", "VAz", "RAz")
    compute_phase_angle(df, "SAz", "VAz", "RAz", "SZA", "VZA", "PA")

    println("Done! Extracted data from ", length(file_list), " files.\n")
    return(df)
end