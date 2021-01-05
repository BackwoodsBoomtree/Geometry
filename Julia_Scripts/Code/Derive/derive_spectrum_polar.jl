using CanopyLayers
using GriddingMachine
using DataFrames

function derive_spectrum_polar(sza::Int64 = 30, clump::Float64 = 1.0)
    FT = Float32;
    ENV["GKSwstype"] = "100";

    # initialize canopy radiation module
    angles, can, can_opt, can_rad, in_rad,
    leaves, rt_con, rt_dim, soil, wls = initialize_rt_module(FT);

    # Initialize output
    sif_682 = zeros(FT, (86, 361));
    sif_757 = zeros(FT, (86, 361));
    sif_771 = zeros(FT, (86, 361));
    ref_682 = zeros(FT, (86, 361));
    ref_757 = zeros(FT, (86, 361));
    ref_771 = zeros(FT, (86, 361));
    rad_682 = zeros(FT, (86, 361));
    rad_757 = zeros(FT, (86, 361));
    rad_771 = zeros(FT, (86, 361));
    ndvi    = zeros(FT, (86, 361));
    nirv    = zeros(FT, (86, 361));
    evi     = zeros(FT, (86, 361));
    lswi    = zeros(FT, (86, 361));

    # Set geometries:
    angles.tts = sza                # SZA

    # Set Clumping
    can.Ω      = clump;

    # Run model
    for psi = 0:360
        angles.psi = psi            # RAA
        for vza = 0:1:85
            angles.tto = vza        # VZA

            # re-run the simulations
            canopy_geometry!(can, angles, can_opt, rt_con);
            canopy_matrices!(leaves, can_opt);
            short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
            canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
            SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

            # Running model with angles of 0, but array rows begin with 1
            sif_682[vza + 1, psi + 1] = can_rad.SIF_obs[8];
            sif_757[vza + 1, psi + 1] = can_rad.SIF_obs[23];
            sif_771[vza + 1, psi + 1] = can_rad.SIF_obs[25] * 0.4667 + can_rad.SIF_obs[26] * 0.5333;

            ref_682[vza + 1, psi + 1] = can_rad.alb_obs[32];
            ref_757[vza + 1, psi + 1] = can_rad.alb_obs[47];
            ref_771[vza + 1, psi + 1] = can_rad.alb_obs[49] * 0.4667 + can_rad.alb_obs[50] * 0.5333;
            
            rad_682[vza + 1, psi + 1] = can_rad.Lo[32];
            rad_757[vza + 1, psi + 1] = can_rad.Lo[47];
            rad_771[vza + 1, psi + 1] = can_rad.Lo[49] * 0.4667 + can_rad.Lo[50] * 0.5333;

            # Vegetation Indices
            nir              = (can_rad.alb_obs[53]  + can_rad.alb_obs[54])  / 2      # 842 and 867 (854.5) in model (MODIS: 841 - 876; 858.5)
            red              = can_rad.alb_obs[25]                                    # 644.5 in model (MODIS: 620 - 670; 645)
            blue             = (can_rad.alb_obs[7]   + can_rad.alb_obs[8])   / 2      # 464.5 and 474.5 (469.5) in model (MODIS: 459 - 479; 469)
            swir             = (can_rad.alb_obs[104] + can_rad.alb_obs[105]) / 2      # 2117 and 2142 (2129.5) in model (MODIS: 2105 - 2155; 2130)
            ndvi[vza + 1, psi + 1]             = (nir - red) ./ (nir + red)
            nirv[vza + 1, psi + 1]             = ndvi[vza + 1, psi + 1] * nir
            evi[vza + 1, psi + 1]              = 2.5 * ((nir - red) / (nir + (6 * red) - (7.5 * blue) + 1))
            lswi[vza + 1, psi + 1]             = (nir - swir) / (nir + swir)
        end
    end

    sif_682_rel = sif_682 ./ rad_682;            # Relative SIF at 682
    sif_757_rel = sif_757 ./ rad_757;            # Relative SIF at 757
    sif_771_rel = sif_771 ./ rad_771;            # Relative SIF at 757
    sif_ratio   = sif_757 ./ sif_771;            # Ratio of SIF 757 and at 771
    sif_nirv    = sif_757 ./ nirv;               # Ratio of SIF 757 and NIRv
    sif_ref     = sif_757 ./ ref_771;

    return sif_682, sif_757, sif_771,
           ref_682, ref_757, ref_771,
           rad_682, rad_757, rad_771,
           sif_682_rel, sif_757_rel, sif_771_rel,
           ndvi, nirv, evi, lswi,
           sif_ref, sif_nirv, sif_ratio,
           wls.WLF, wls.WL
end