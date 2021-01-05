using Plots
using CanopyLayers
using GriddingMachine
using DataFrames
using CSV

function derive_spectrum_psi(psi::Int64 = 180, LAI = 4, clump = 1.0)
    FT = Float32;
    ENV["GKSwstype"] = "100";

    # initialize canopy radiation module
    angles, can, can_opt, can_rad, in_rad,
	leaves, rt_con, rt_dim, soil, wls = initialize_rt_module(FT);

	# Initialize output
	sif_757_nadir = FT[];
	sif_757_glint = FT[];
	sif_hemi      = FT[];

	angles.psi 		= psi            # RAA
	can.LAI 		= LAI
	can.Ω           = clump;
	println("RAA = ", RAA, "; LAI = ", LAI, "; Clumping Index = ", clump)

	# Run model
	for tts = 0:75              # SZA
		
		# Run for NADIR where VZA = 0
		angles.tts = tts        # SZA
		angles.tto = 0          # VZA

		# re-run the simulations
		canopy_geometry!(can, angles, can_opt, rt_con);
		canopy_matrices!(leaves, can_opt);
		short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
		canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
		SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

		push!(sif_757_nadir, can_rad.SIF_obs[23]);
		push!(sif_hemi,      can_rad.SIF_hemi[23]);

		# Run for Glint
		angles.tto = tts - 10*sind(tts)          # VZA for Glint

		# re-run the simulations
		canopy_geometry!(can, angles, can_opt, rt_con);
		canopy_matrices!(leaves, can_opt);
		short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
		canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves, wls, rt_con);
		SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);

		push!(sif_757_glint, can_rad.SIF_obs[23]);
	end
	return sif_757_nadir, sif_757_glint, sif_hemi
end

RAA   = 180
LAI   = 4
clump = 0.4

sif_757_nadir, sif_757_glint, sif_hemi = derive_spectrum_psi(RAA, LAI, clump);

sif_avg_total    = sif_hemi ./ π; # Average Outgoing SIF
nadir_avg_total  = sif_757_nadir ./ sif_avg_total;
glint_avg_total  = sif_757_glint ./ sif_avg_total;
nadir_glint      = sif_757_nadir ./ sif_757_glint;

plot(sif_757_nadir, label = "NADIR (VZA = 0°)", legend = :topleft)
plot!(sif_757_glint, label = "Glint    (VZA = SZA - 10*sind(SZA))")
plot!(sif_avg_total, label = "Total Average Outgoing SIF (mW m⁻² nm⁻¹)")
plot!(nadir_avg_total, label = "NADIR / Total")
plot!(glint_avg_total, label = "Glint / Total")
plot!(nadir_glint, label = "NADIR / Glint")
xlabel!("Solar Zenith Angle (°)")
ylabel!("SIF₇₅₇ (mW m⁻² nm⁻¹ sr⁻¹)")
title!(string("RAA = ", RAA, ", LAI = ", LAI, ", CI = ", clump))

sif_757_nadir_norm = (sif_757_nadir .- minimum(sif_757_nadir)) ./ (maximum(sif_757_nadir) - minimum(sif_757_nadir));
sif_757_glint_norm = (sif_757_glint .- minimum(sif_757_glint)) ./ (maximum(sif_757_glint) - minimum(sif_757_glint));
sif_avg_total_norm = (sif_avg_total .- minimum(sif_avg_total)) ./ (maximum(sif_avg_total) - minimum(sif_avg_total));

nadir_avg_total_norm = (nadir_avg_total .- minimum(nadir_avg_total)) ./ (maximum(nadir_avg_total) - minimum(nadir_avg_total));
glint_avg_total_norm = (glint_avg_total .- minimum(glint_avg_total)) ./ (maximum(glint_avg_total) - minimum(glint_avg_total));
nadir_glint_norm     = (nadir_glint .- minimum(nadir_glint)) ./ (maximum(nadir_glint) - minimum(nadir_glint));

plot(sif_757_nadir_norm, label = "NADIR (VZA = 0°)", legend = :topleft)
plot!(sif_757_glint_norm, label = "Glint    (VZA = SZA - 10*sind(SZA))")
plot!(sif_avg_total_norm, label = "Total Average Outgoing SIF (mW m⁻² nm⁻¹)")
plot!(nadir_avg_total_norm, label = "NADIR / Total")
plot!(glint_avg_total_norm, label = "Glint / Total")
plot!(nadir_glint_norm, label = "NADIR / Glint")
xlabel!("Solar Zenith Angle (°)")
ylabel!("SIF₇₅₇ (mW m⁻² nm⁻¹ sr⁻¹)")
title!(string("RAA = ", RAA, ", LAI = ", LAI, ", CI = ", clump))

df = DataFrame(
	sif_757_nadir   = sif_757_nadir,
	sif_757_glint   = sif_757_glint,
	sif_avg_total   = sif_avg_total,
	nadir_avg_total = nadir_avg_total,
	glint_avg_total = glint_avg_total,
	nadir_glint     = nadir_glint,
	sif_757_nadir_norm   = sif_757_nadir_norm,
	sif_757_glint_norm   = sif_757_glint_norm,
	sif_avg_total_norm   = sif_avg_total_norm,
	nadir_avg_total_norm = nadir_avg_total_norm,
	glint_avg_total_norm = glint_avg_total_norm,
	nadir_glint_norm     = nadir_glint_norm
	)

CSV.write("C:/Russell/Projects/Geometry/Julia_Scripts/CSV/proposal_plots_BRDF_Clumping_0.4.csv", df)