
### Computing Relative Azimuth and Phase Angles ###

using DataFrames

function compute_relative_angle(df::DataFrame, saz_name::String, vaz_name::String, raz_name::String) # Column names for input, and output raz
    # Calculate Relative Azimuth Angle
    df[Symbol(raz_name)] = broadcast(abs, (df[Symbol(saz_name)] - df[Symbol(vaz_name)]))
    for a in 1:length(df[Symbol(raz_name)])
        if df[Symbol(raz_name)][a] > 180
            df[Symbol(raz_name)][a] = abs(df[Symbol(raz_name)][a] - 360)
        end
    end
    println("Computed relative azimuth angles.")
    return(df)
end

function compute_phase_angle(df::DataFrame, saz_name::String, vaz_name::String, raz_name::String, sza_name::String, vza_name::String, pa_name::String) # Column names for input, and output pa
    ## Phase angles are set to negative values if the observational azimuth angle is bigger than the solar azimuth angle
    ## (negative phase angle means the sun is to the right of the satellite)
    df[Symbol(pa_name)] = 0.0
    for i in 1:length(df[Symbol(saz_name)])
        if df[Symbol(vaz_name)][i] > df[Symbol(saz_name)][i]
            df[Symbol(pa_name)][i] = -1
        elseif df[Symbol(vaz_name)][i] < df[Symbol(saz_name)][i]
            df[Symbol(pa_name)][i] = 1
        end
    end
    for i in 1:length(df[Symbol(pa_name)])
        df[Symbol(pa_name)][i] = (acos(cosd(df[Symbol(sza_name)][i]) * cosd(df[Symbol(vza_name)][i]) + sind(df[Symbol(vza_name)][i]) * sind(df[Symbol(sza_name)][i]) * cosd(df[Symbol(raz_name)][i])) * 180 / π) * df[Symbol(pa_name)][i]
    end
    println("Computed phase angles.")
    return(df)
end