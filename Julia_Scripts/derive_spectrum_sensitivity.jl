using CanopyLayers
using GriddingMachine
using DataFrames

struct SIFComparison2020{FT} end

function derive_spectrum_sensitivity(model_variable::String, var_range::Vector{}, by)
    FT = Float32;
    ENV["GKSwstype"]="100";

    # initialize canopy radiation module
    angles, can, can_opt, can_rad, in_rad,
    leaves, rt_con, rt_dim, soil, wls = initialize_rt_module(FT);

    # create a copy of incoming radiation
    in_rad_bak = deepcopy(in_rad);
    e_all_dire = sum(in_rad_bak.E_direct  .* wls.dWL) / 1000;
    e_all_diff = sum(in_rad_bak.E_diffuse .* wls.dWL) / 1000;

    # create a matrix to store the spectrum
    var_values = collect(var_range[1]:by:var_range[2])
    mat_REF = zeros(FT, (length(var_values), length(wls.WL)));
    mat_SIF = zeros(FT, (length(var_values), length(wls.WLF)));

    # # Load LUTs
    # if source_cab == "LUT"
    #     CAB_LUT = load_LUT(LeafChlorophyll{FT}());
    #     mask_LUT!(CAB_LUT, FT[0,Inf]);
    # end

    # iterate through the data
    for i in 1:length(var_values)
        if model_variable == "Ω"
            println("Canopy Clumping Index set to ", var_values[i])
            can.Ω = var_values[i];
        elseif model_variable == "LAI"
            can.LAI = var_values[i];
            can.iLAI = can.LAI * can.Ω / can.nLayer;
        # elseif var == "rad"
        #     in_rad = deepcopy(in_rad_bak);
        #     in_rad.E_direct  .*= input_data.incoming_direct_era5[i]  / e_all_dire;
        #     in_rad.E_diffuse .*= input_data.incoming_diffuse_era5[i] / e_all_diff;
        elseif model_variable == "tts" # SZA
            angles.tts = var_values[i]
        elseif model_variable == "psi" # RAz
            angles.psi = var_values[i]
        elseif model_variable == "tto" # VZA
            angles.tto = var_values[i]
        end

        # Run fluspect on each layer to change Cab
        for j in 1:20
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
    output_data.temp = var_values
    names!(output_data, Symbol.([model_variable]))
    output_data.SIF740 = mat_SIF[:,19] .* 0.4 .+ mat_SIF[:,20] .* 0.6;
    output_data.SIF757 = mat_SIF[:,23];
    output_data.SIF771 = mat_SIF[:,25] .* 0.4667 .+ mat_SIF[:,26] .* 0.5333;
    output_data.REF757 = mat_REF[:,47];
    output_data.REF771 = mat_REF[:,49] .* 0.4667 .+ mat_REF[:,50] .* 0.5333;
    output_data.SIF757_Relative = output_data.SIF757 ./ output_data.REF757
    output_data.SIF771_Relative = output_data.SIF771 ./ output_data.REF771

    return mat_REF, mat_SIF, wls.WLF, wls.WL, output_data
end