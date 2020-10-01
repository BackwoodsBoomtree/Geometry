#
# This function requires the following packages
#
# CSV
# DataFrames
# CanopyLayers
#
using CSV
using CanopyLayers
using DataFrames
using PlotPlants

struct SIFComparison2020{FT} end

"""
    derive_spectrum()

Derive and plot the albedo and SIF spectrum
"""
function derive_spectrum()
    FT = Float32;

    # read data from file
    # TODO add clumping factor
    # TODO add cab
    _file = joinpath(@__DIR__, "../../data/russ_data.csv");
    _data = DataFrame!(CSV.File(_file));

    # initialize canopy radiation module
    angles, can, can_opt, can_rad, in_rad,
    leaves, rt_con, rt_dim, soil, wls = initialize_rt_module(FT);

    # create a copy of incoming radiation
    in_rad_bak = deepcopy(in_rad);
    e_all_dire = sum(in_rad_bak.E_direct  .* wls.dWL) / 1000;
    e_all_diff = sum(in_rad_bak.E_diffuse .* wls.dWL) / 1000;

    # create a matrix to store the spectrum
    mat_REF = zeros(FT, (length(_data.vza), length(wls.WL )));
    mat_SIF = zeros(FT, (length(_data.vza), length(wls.WLF)));

    # iterate through the data
    for i in eachindex(_data.vza)
        # change the angles
        angles.tts = _data.sza[i];
        angles.psi = _data.raa[i];
        angles.tto = _data.vza[i];

        # change the canopy profiles
        can.LAI    = _data.lai_era5[i];
        can.iLAI   = _data.lai_era5[i] * can.Ω / can.nLayer;

        # change incoming radiation, comment these to see the improvements
        in_rad = deepcopy(in_rad_bak);
        in_rad.E_direct  .*= _data.incoming_direct_era5[i]  / e_all_dire;
        in_rad.E_diffuse .*= _data.incoming_diffuse_era5[i] / e_all_diff;

        # re-run the simulations
        canopy_geometry!(can, angles, can_opt, rt_con);
        canopy_matrices!(leaves, can_opt);
        short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
        canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
        SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

        mat_REF[i,:] .= can_rad.alb_obs;
        mat_SIF[i,:] .= can_rad.SIF_obs;
    end

    _fig1,_axes = create_canvas(1; figsize=(7.5,3.5));
    _ax1 = _axes[1];
    _cm1 = _ax1.pcolor(mat_REF);
    _ax1.set_title("Reflectance Spetrum");
    _ax1.set_xlabel("Wave Length Index");
    _ax1.set_ylabel("Site Index");
    _fig1.colorbar(_cm1, ax=_ax1, label="Albedo");
    _fig1.set_tight_layout(true);

    _fig2,_axes = create_canvas(2; figsize=(5.5,3.5));
    _ax2 = _axes[1];
    _cm2 = _ax2.pcolor(wls.WLF, eachindex(_data.vza), mat_SIF);
    _ax2.set_title("SIF Spetrum");
    _ax2.set_xlabel("Wave Length (nm)");
    _ax2.set_ylabel("Site Index");
    _fig2.colorbar(_cm2, ax=_ax2, label="SIF");
    _fig2.set_tight_layout(true);

    # 742 nm and 737 nm to 740 nm
    # 757 nm to 757 nm
    # 767 nm and 774.5 nm to 771 nm
    _sif740 = mat_SIF[:,19] .* 0.4 .+ mat_SIF[:,20] .* 0.6;
    _sif757 = mat_SIF[:,23];
    _sif771 = mat_SIF[:,25] .* 0.4667 .+ mat_SIF[:,26] .* 0.5333;

    _fig3,_axes = create_canvas(3; figsize=(3.5,3.5));
    _ax3 = _axes[1];
    _ax3.plot(_sif740, _data.sif740, "k+");
    plot_line_regress(_ax3, _sif740, _data.sif740; interval=true);
    _ax3.set_xlabel("SIF 740");
    _ax3.set_ylabel("SIF 740");
    _fig3.set_tight_layout(true);

    _fig4,_axes = create_canvas(4; figsize=(3.5,3.5));
    _ax4 = _axes[1];
    _ax4.plot(_sif757, _data.sif757, "k+");
    plot_line_regress(_ax4, _sif757, _data.sif757; interval=true);
    _ax4.set_xlabel("SIF 757");
    _ax4.set_ylabel("SIF 757");
    _fig4.set_tight_layout(true);

    _fig5,_axes = create_canvas(5; figsize=(3.5,3.5));
    _ax5 = _axes[1];
    _ax5.plot(_sif771, _data.sif771, "k+");
    plot_line_regress(_ax5, _sif771, _data.sif771; interval=true);
    _ax5.set_xlabel("SIF 771");
    _ax5.set_ylabel("SIF 771");
    _fig5.set_tight_layout(true);

    return nothing
end
