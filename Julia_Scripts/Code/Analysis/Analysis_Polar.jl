code_dir = "C:/Russell/Projects/Geometry/Julia_Scripts/Code/"
include(code_dir*"Plot/Plot_Polar.jl")

#### Regular Plot ####
sza       = 30;
vza_max   = 85;
plot_type = "757_clump"                  # plot_type is "SRRR" or "757_vis" or "757_clump"
clump     = 0.4
file_name = "C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Polar/proposal_BRDF.pdf"

p = plot_polar(sza, vza_max, plot_type, clump)
p.savefig(file_name)

#### Animated Plot ####
dir_temp_anim  = "C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Polar/animate_temp/"
file_name_anim = "C:/Russell/Projects/Geometry/Julia_Scripts/Figures/Polar/animated_BRDF.gif"
plot_polar_anim(dir_temp_anim, file_name_anim)