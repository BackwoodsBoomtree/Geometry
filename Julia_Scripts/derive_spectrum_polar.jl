using CanopyLayers
using GriddingMachine
using DataFrames

struct SIFComparison2020{FT} end

function derive_spectrum(sza)
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

    # Test a VZA dependence in the principal plane ####
    SIF_FR = Float32[]
    SIF_R = Float32[]
    reflVIS = Float32[]
    reflNIR = Float32[]

    # Set geometries:
    angles.tts = sza                # SZA
    angles.psi = 0                  # Relative azimuth angle (principal plane)
    VZA = collect(-89.5:0.5:89.5)   # VZA range

    for psi = 0:360
        angles.psi = psi
        for VZA = 0:1:85
            angles.tto = VZA
    
            compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
            compute_canopy_matrices!(arrayOfLeaves, canOpt_rt);
            simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil);
            computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set);

            # re-run the simulations
            canopy_geometry!(can, angles, can_opt, rt_con);
            canopy_matrices!(leaves, can_opt);
            short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
            canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
            SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

            push!(SIF_FR, can_rad.SIF_obs[20])
            push!(SIF_R , can_rad.SIF_obs[8])
            push!(reflVIS, can_rad.alb_obs[28])
            push!(reflNIR, can_rad.alb_obs[52])
        end
    end

    return SIF_FR, SIF_R, reflVIS, reflNIR, wls.WLF, wls.WL
end