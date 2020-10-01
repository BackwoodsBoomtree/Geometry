#
# This function requires the following packages
#
# CSV
# DataFrames
# CanopyLayers
# PlotPlants (optional, if you do not use my plotting commands)
#

struct SIFComparison2020{FT} end

"""
    derive_spectrum(proj::SIFComparison2020{FT};
                    saving::Bool,
                    use_latex::Bool) where {FT<:AbstractFloat}

Derive and plot the albedo and SIF spectrum, given
- `proj` [`SIFComparison2020`](@ref) type project control
- `saving` Optional. If true, save the figure
- `use_latex` Optional. If true, use latex and serif font
"""
function derive_spectrum(
            proj::SIFComparison2020{FT};
            saving::Bool = false,
            use_latex::Bool = true
) where {FT<:AbstractFloat}
    # use latex and serif font
    if use_latex use_serif_tex(); end

    # read data from file
    # TODO add clumping factor
    # TODO add cab
    _file = joinpath(@__DIR__, "../../data/russ_data.csv");
    _data = DataFrame!(CSV.File(_file));

    # initialize canopy radiation module
    angles, can, can_opt, can_rad, in_rad,
    leaves, rt_con, rt_dim, soil, wls = initialize_rt_module(FT);

    # create a matrix to store the spectrum
    mat_REF = zeros(FT, (length(_data.vza), length(wls.WL )));
    mat_SIF = zeros(FT, (length(_data.vza), length(wls.WLF)));

    # iterate through the data
    for i in eachindex(_data.vza)
        angles.tts = _data.sza[i];
        angles.psi = _data.raa[i];
        angles.tto = _data.vza[i];
        can.LAI    = _data.lai_era5[i];
        can.iLAI   = _data.lai_era5[i] * can.Ω / can.nLayer;

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

    if saving
        _fig1.savefig("figures/2020_sif_comparison/ref_spectrum.pdf",
                      bbox_inches="tight");
        _fig2.savefig("figures/2020_sif_comparison/sif_spectrum.pdf",
                      bbox_inches="tight");
    end

    return nothing
end
