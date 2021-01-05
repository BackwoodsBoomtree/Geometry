include("derive_spectrum_polar.jl")
using Plots
using FileIO, ImageMagick, Colors, FixedPointNumbers

x = collect(-180:180)
y = collect(0:85)
dir_temp_anim  = "C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Polar/animate_temp/"
file_name_anim = "C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Polar/animated_3D_BRDF.gif"

pyplot()
n = 1;
for sza = 0:1:85
	# Create layer
	z = fill(sza, (86, 361));

	# Run model
	sif_682, sif_757, sif_771,
	ref_682, ref_757, ref_771,
	rad_682, rad_757, rad_771,
	sif_682_rel, sif_757_rel, sif_771_rel,
	ndvi, nirv, evi, lswi,
	sif_ref, sif_nirv, sif_ratio,
	wlf, wl = derive_spectrum_polar(sza);

	# Reorder SIF array for -180 to 180 RAA
	sif_757_r = hcat(sif_757[:, 181:361], sif_757[:, 1:180])

	# Plot
	if sza == 0
		# Make Dir
		if isdir(dir_temp_anim) == false # Create directory
        	mkdir(dir_temp_anim)
    	else
			rm(dir_temp_anim, recursive = true) # if exists, delete and recreate to avoid duplicating images
			mkdir(dir_temp_anim)
		# Plot
		p = surface(x, y, z, zlim = (0, 85), fill_z = sif_757_r, colorbar = true, colorbar_title = "\n", surfacealpha = 1,
					ylabel = "Viewing Zenith Angle (°)", xlabel = "Relative Azimuth Angle (°)",
					zlabel = "Solar Zenith Angle (°)", zguidefontrotation = 90,
					yguidefontrotation = -46, xguidefontrotation = 18, clims = (0, 2.5),
					title = "SIF₇₅₇", size = (700,600), camera = (-30,30))
		end
	else
		p = surface!(x,y,z, fill_z = sif_757_r, colorbar = true)
	end
	savefig(p, dir_temp_anim*"temp"*string(lpad(n, 3, '0'))*".png")
	n = n + 1
end

# Make animation
# Create an image stack for animation
n_img = length(readdir(dir_temp_anim; join = true));
img = ()
for i in 1:n_img
	img_path = dir_temp_anim*"temp"*string(lpad(i, 3, '0'))*".png"
	if i == 1
		img = load(img_path)
	else
		img = cat(img, load(img_path), dims=3)
	end
end
# Save image stack as animated GIF
save(file_name_anim, img; fps = 15)
# Remove temp directory
rm(dir_temp_anim, recursive = true)