using CanopyLayers
using GriddingMachine
using DataFrames

struct SIFComparison2020{FT} end

function derive_spectrum(input_data, source_cab, week, clumping)
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
    mat_REF = zeros(FT, (length(input_data.VZA), length(wls.WL)));
    mat_SIF = zeros(FT, (length(input_data.VZA), length(wls.WLF)));

    # Load LUTs
    if source_cab == "LUT"
        CAB_LUT = load_LUT(LeafChlorophyll{FT}());
        mask_LUT!(CAB_LUT, FT[0,Inf]);
    end

    # iterate through the data
    for i in eachindex(input_data.VZA)
        # change the angles
        angles.tts = input_data.SZA[i];
        angles.psi = input_data.RAz[i];
        angles.tto = input_data.VZA[i];

        # change the canopy profiles
        if clumping == "clump"
            can.Ω = input_data.clump[i];
        else
            can.Ω = 1
        end
        can.LAI = input_data.lai_cop[i];
        can.iLAI = can.LAI * can.Ω / can.nLayer;

        # change incoming radiation
        in_rad = deepcopy(in_rad_bak);
        in_rad.E_direct  .*= input_data.incoming_direct_era5[i]  / e_all_dire;
        in_rad.E_diffuse .*= input_data.incoming_diffuse_era5[i] / e_all_diff;

        # Run fluspect on each layer to change Cab
        for j in 1:20
            # SIFyield
            leaves[j].fqe = 0.47
            leaves[j].Cx = 0.22

            if source_cab == "LUT"
                leaves[j].Cab = read_LUT(CAB_LUT, FT(input_data.latitude[i]), FT(input_data.longitude[i]), week);
            else
                leaves[j].Cab = source_cab;
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
    output_data = deepcopy(input_data)
    output_data.SIF740_Sim = mat_SIF[:,19] .* 0.4 .+ mat_SIF[:,20] .* 0.6;
    output_data.SIF757_Sim = mat_SIF[:,23];
    output_data.SIF771_Sim = mat_SIF[:,25] .* 0.4667 .+ mat_SIF[:,26] .* 0.5333;
    output_data.REF757_Sim = mat_REF[:,47];
    output_data.REF771_Sim = mat_REF[:,49] .* 0.4667 .+ mat_REF[:,50] .* 0.5333;

    return mat_REF, mat_SIF, wls.WLF, wls.WL, output_data
end