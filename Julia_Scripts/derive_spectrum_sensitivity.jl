using CanopyLayers
using GriddingMachine
using DataFrames

struct SIFComparison2020{FT} end

function derive_spectrum_sensitivity(model_variable::String, var_range::Array{Float64,1}, by)
    FT = Float32;
    ENV["GKSwstype"]="100";

    # create vector of values within input range and by input step
    var_values = collect(var_range[1]:by:var_range[2])
    
    if model_variable == "rad_ratio"
        var_values_front = deepcopy(var_values)
        var_values_front = reverse(deleteat!(var_values_front,1))
        var_values_front_rat = inv.(rationalize.(var_values_front))
        var_values_rat = rationalize.(var_values)
        var_values_front = round.(inv.(var_values_front), digits = 2)
        var_values = [var_values_front; var_values]
        var_values_rat = [var_values_front_rat; var_values_rat]
    end

    # initialize canopy radiation module
    angles, can, can_opt, can_rad, in_rad,
    leaves, rt_con, rt_dim, soil, wls = initialize_rt_module(FT);

    # create a copy of incoming radiation
    in_rad_bak = deepcopy(in_rad);
    e_all_dire = sum(in_rad_bak.E_direct  .* wls.dWL) / 1000;
    e_all_diff = sum(in_rad_bak.E_diffuse .* wls.dWL) / 1000;

    # create a matrix to store the spectrum
    mat_REF = zeros(FT, (length(var_values), length(wls.WL)));
    mat_SIF = zeros(FT, (length(var_values), length(wls.WLF)));

    # Run the model for each step in the range
    for i in 1:length(var_values)
        
        ### Canopy Variables ###
        if model_variable == "nLayer"           # Canopy Layers (20)
            angles, can, can_opt, can_rad, in_rad,
            leaves, rt_con, rt_dim, soil, wls = initialize_rt_module(FT, nLayer = var_values[i]);
        elseif model_variable == "Ω"            # Clumping Index (1.0)
            can.Ω             = var_values[i];
            can.iLAI          = can.LAI * can.Ω / can.nLayer;
        elseif model_variable == "LAI"          # Leaf Area Index (3.0)
            can.LAI           = var_values[i];
            can.iLAI          = can.LAI * can.Ω / can.nLayer;
        elseif model_variable == "clump_a"      # Structure Factor a (1.0)
            can.clump_a       = var_values[i];
        elseif model_variable == "clump_b"      # Structure Factor b (0.0)
            can.clump_b       = var_values[i];
        elseif model_variable == "leaf_width"   # Leaf Width (0.1)
            can.leaf_width    = var_values[i];
        elseif model_variable == "hc"           # Vegetation Height (2.0)
            can.hc            = var_values[i];
        elseif model_variable == "LIDFa"        # Leaf Inclination (0.0)
            can.LIDFa         = var_values[i];
        elseif model_variable == "LIDFb"        # Variation in leaf Inclination (0.0)
            can.LIDFb         = var_values[i];
        elseif model_variable == "hot"          # HotSpot Parameter (still need to check!) (0.05)
            can.hot           = var_values[i];
        elseif model_variable == "height"       # Canopy height (20 m)
            can.height        = var_values[i];
        elseif model_variable == "z0m"          # Canopy roughness (1 m)
            can.z0m           = var_values[i];
        elseif model_variable == "z0h"          # Tree roughness (-999.0 m)
            can.z0h           = var_values[i];
        elseif model_variable == "d"            # Canopy displacement height (-999.0 m)
            can.d             = var_values[i];
        elseif model_variable == "Cd"           # Turbulent Transfer Coefficient m/sqrt(s) (0.01)
            can.Cd            = var_values[i];
        end

        ### Incoming Radiation ###
        if model_variable == "rad_ratio"        # Ratio of direct to diffuse light Radiation
            in_rad             = deepcopy(in_rad_bak);
            direct_incoming    = (var_values[i] / (var_values[i] + 1)) * 500
            diffuse_incoming   = 500 - direct_incoming
            in_rad.E_direct  .*= direct_incoming  / e_all_dire;
            in_rad.E_diffuse .*= diffuse_incoming / e_all_diff;
            println("direct = ", round(direct_incoming, digits = 3),
                    " and diffuse = ", round(diffuse_incoming, digits = 3))
        end

        ### Sun Sensor Geometry Variables ###
        if model_variable == "tts"              # Solar Zenith Angle (30 degrees)
            angles.tts    = var_values[i]
        elseif model_variable == "tto"          # Viewing Zenith Angle (0 degrees)
            angles.tto        = var_values[i]
        elseif model_variable == "psi"          # Relative Azimuth Angle (0 degrees)
            angles.psi        = var_values[i]
        end

        println(model_variable, " set to ", var_values[i])

        # Run fluspect on each layer
        for j in 1:can.nLayer
            if model_variable == "fqe" # SIFyield
                leaves[j].fqe = var_values[i]
            elseif model_variable == "Cx"
                leaves[j].Cx = var_values[i]
            elseif model_variable == "Cab" # Chlorophyll ab
                leaves[j].Cab = var_values[i]
            end
            fluspect!(leaves[j], wls);
        end

        # re-run the simulations
        canopy_geometry!(can, angles, can_opt, rt_con);
        canopy_matrices!(leaves, can_opt);
        short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
        canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
        SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

        mat_REF[i,:] .= can_rad.Lo;
        mat_SIF[i,:] .= can_rad.SIF_obs;
    end

    # 742 nm and 737 nm to 740 nm
    # 757 nm to 757 nm
    # 767 nm and 774.5 nm to 771 nm
    output_data = DataFrame()
    if model_variable == "rad_ratio"
        output_data.temp = var_values_rat
    else
        output_data.temp = var_values
    end
    names!(output_data, Symbol.([model_variable]));
    output_data.SIF740 = mat_SIF[:,19] .* 0.4 .+ mat_SIF[:,20] .* 0.6;
    output_data.SIF757 = mat_SIF[:,23];
    output_data.SIF771 = mat_SIF[:,25] .* 0.4667 .+ mat_SIF[:,26] .* 0.5333;
    output_data.REF757 = mat_REF[:,47];
    output_data.REF771 = mat_REF[:,49] .* 0.4667 .+ mat_REF[:,50] .* 0.5333;
    output_data.SIF757_Relative = output_data.SIF757 ./ output_data.REF757
    output_data.SIF771_Relative = output_data.SIF771 ./ output_data.REF771

    return mat_REF, mat_SIF, wls.WLF, wls.WL, output_data
end