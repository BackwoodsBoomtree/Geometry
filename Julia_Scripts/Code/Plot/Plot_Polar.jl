code_dir = "C:/Russell/Projects/Geometry/Julia_Scripts/Code/"
include(code_dir*"Derive/derive_spectrum_polar.jl")
using PyPlot
using FileIO, ImageMagick, Colors, FixedPointNumbers

function plot_polar(sza, vza_max, plot_type, clump)
    sif_682, sif_757, sif_771,
    ref_682, ref_757, ref_771,
    rad_682, rad_757, rad_771,
    sif_682_rel, sif_757_rel, sif_771_rel,
    ndvi, nirv, evi, lswi,
    sif_ref, sif_nirv, sif_ratio,
    wlf, wl = derive_spectrum_polar(sza, 1.0);
    
    if plot_type == "SRRR"              # SIF, Radiance, Relative SIF, Reflectance at 682, 757, 771
                                        # Set up main plot for SIF, Radiance, SIF rel, Reflectance
        p = figure(figsize = (14, 16))  # w, h
        suptitle(string("BRDF for SZA = ", sza), fontsize = 16, y = 0.935)
        p.subplots_adjust(hspace = 0.4) # horizontal spacing between subplots
                                            # SIF PLOTS
                                            # SIF 682
        subplot(4,3,1, polar = true)
        grid(false)
        plot_sif_682 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_682[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₆₈₂", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11) 
                                            # SIF 757
        subplot(4,3,2, polar = true)
        grid(false)
        plot_sif_757 = PyPlot.contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # SIF 771
        subplot(4,3,3, polar = true)
        grid(false)
        plot_sif_771 = PyPlot.contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_771[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₇₁", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Radiance PLOTS
                                            # Radiance 682
        subplot(4,3,4, polar = true)
        grid(false)
        plot_rad_682 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), rad_682[1:(vza_max + 1), :], cmap = :viridis)
        title("Radiance₆₈₂", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Radiance 757
        subplot(4,3,5, polar = true)
        grid(false)
        plot_rad_757 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), rad_757[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Radiance₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Radiance 771
        subplot(4,3,6, polar = true)
        grid(false)
        plot_rad_771 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), rad_771[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Radiance₇₇₁", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Relative SIF PLOTS
                                            # Relative SIF 682
        subplot(4,3,7, polar = true)
        grid(false)
        plot_sif_682_rel = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_682_rel[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Relative SIF₆₈₂", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Relative SIF 757
        subplot(4,3,8, polar = true)
        grid(false)
        plot_sif_757_rel = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757_rel[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Relative SIF₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Relative SIF 771
        subplot(4,3,9, polar = true)
        grid(false)
        plot_sif_771_rel = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_771_rel[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Relative SIF₇₇₁", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Reflectance PLOTS
                                            # Reflectance 682
        subplot(4,3,10, polar = true)
        grid(false)
        plot_ref_682 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), ref_682[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Reflectance₆₈₂", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Reflectance 771
        subplot(4,3,12, polar = true)
        grid(false)
        plot_ref_771 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), ref_771[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Reflectance₇₇₁", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)

        close(p)

    elseif plot_type == "757_vis"
                                            # Plot with VIs
        p = figure(figsize = (14, 13))  # w, h
        suptitle(string("BRDF for SZA = ", sza), fontsize = 16, y = 0.935)
        p.subplots_adjust(hspace = 0.4) # horizontal spacing between subplots
                                            # SIF 757
        subplot(3,3,1, polar = true)
        grid(false)
        plot_sif_757 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Radiance 757
        subplot(3,3,2, polar = true)
        grid(false)
        plot_rad_757 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), rad_757[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Radiance₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Reflectance 757
        subplot(3,3,3, polar = true)
        grid(false)
        plot_ref_757 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), ref_757[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Reflectance₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # SIF/NIRv
        subplot(3,3,4, polar = true)
        grid(false)
        plot_sif_nirv = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_nirv[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₅₇/NIRv", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Relative SIF 757
        subplot(3,3,5, polar = true)
        grid(false)
        plot_sif_757_rel = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757_rel[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Relative SIF₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # SIF/Ref 757
        subplot(3,3,6, polar = true)
        grid(false)
        plot_sif_ref = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_ref[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₅₇/Reflectance₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # NIRv
        subplot(3,3,7, polar = true)
        grid(false)
        plot_nirv = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), nirv[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("NIRv", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # EVI
        subplot(3,3,8, polar = true)
        grid(false)
        plot_evi = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), evi[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("EVI", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # NDVI
        subplot(3,3,9, polar = true)
        grid(false)
        plot_ndvi = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), ndvi[1:(vza_max + 1), :], cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("NDVI", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)

        close(p)

    elseif plot_type == "757_vis_animated"  # Animated Plot with VIs
                                            
        p = figure(figsize = (14, 13))  # w, h
        suptitle(string("BRDF for SZA = ", sza), fontsize = 16, y = 0.935)
        p.subplots_adjust(hspace = 0.4) # horizontal spacing between subplots
                                            # SIF 757
        subplot(3,3,1, polar = true)
        grid(false)
        plot_sif_757 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757[1:(vza_max + 1), :], cmap = :viridis, levels = collect(1.86:0.06:2.34), extend = "both")
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Radiance 757
        subplot(3,3,2, polar = true)
        grid(false)
        plot_rad_757 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), rad_757[1:(vza_max + 1), :], cmap = :viridis, levels = collect(97.5:2.5:120), extend = "both")
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Radiance₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Reflectance 757
        subplot(3,3,3, polar = true)
        grid(false)
        plot_ref_757 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), ref_757[1:(vza_max + 1), :], cmap = :viridis, levels = collect(0.37:0.01:0.45), extend = "both")
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Reflectance₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # SIF/NIRv
        subplot(3,3,4, polar = true)
        grid(false)
        plot_sif_nirv = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_nirv[1:(vza_max + 1), :], cmap = :viridis, levels = collect(5.25:0.05:5.70), extend = "both")
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₅₇/NIRv", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Relative SIF 757
        subplot(3,3,5, polar = true)
        grid(false)
        plot_sif_757_rel = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757_rel[1:(vza_max + 1), :], cmap = :viridis, levels = collect(0.01875:0.00015:0.02010), extend = "both")
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Relative SIF₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # SIF/Ref 757
        subplot(3,3,6, polar = true)
        grid(false)
        plot_sif_ref = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_ref[1:(vza_max + 1), :], cmap = :viridis, levels = collect(4.70:0.05:5.10), extend = "both")
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₅₇/Reflectance₇₅₇", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # NIRv
        subplot(3,3,7, polar = true)
        grid(false)
        plot_nirv = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), nirv[1:(vza_max + 1), :], cmap = :viridis, levels = collect(0.352:0.008:0.424), extend = "both")
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("NIRv", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # EVI
        subplot(3,3,8, polar = true)
        grid(false)
        plot_evi = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), evi[1:(vza_max + 1), :], cmap = :viridis, levels = collect(0.67:0.01:0.76), extend = "both")
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("EVI", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # NDVI
        subplot(3,3,9, polar = true)
        grid(false)
        plot_ndvi = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), ndvi[1:(vza_max + 1), :], cmap = :viridis, levels = collect(0.876:0.004:0.908), extend = "both")
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("NDVI", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)

        close(p)
    elseif plot_type == "757_clump"
        sif_682_c, sif_757_c, sif_771_c,
        ref_682_c, ref_757_c, ref_771_c,
        rad_682_c, rad_757_c, rad_771_c,
        sif_682_rel_c, sif_757_rel_c, sif_771_rel_c,
        ndvi_c, nirv_c, evi_c, lswi_c,
        sif_ref_c, sif_nirv_c, sif_ratio_c,
        wlf, wl = derive_spectrum_polar(sza, clump);

        p = figure(figsize = (14, 13))  # w, h
        suptitle(string("BRDF for SZA = ", sza), fontsize = 16, y = 0.935)
        p.subplots_adjust(hspace = 0.4) # horizontal spacing between subplots
                                            # SIF 757
        subplot(3,3,1, polar = true)
        grid(false)
        plot_sif_757 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757[1:(vza_max + 1), :], levels = collect(1.2:0.2:2.6), extend = "both", cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₅₇ Ω = 1", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # NIRv
        subplot(3,3,2, polar = true)
        grid(false)
        plot_nirv = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), nirv[1:(vza_max + 1), :], levels = collect(0.20:0.03:0.44), extend = "both", cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("NIRv Ω = 1", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Relative SIF 757
        subplot(3,3,3, polar = true)
        grid(false)
        plot_sif_757_rel = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757_rel[1:(vza_max + 1), :], levels = collect(0.0172:0.0004:0.0208), extend = "both", cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Relative SIF₇₅₇ Ω = 1", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # SIF 757 Ω = 0.4
        subplot(3,3,4, polar = true)
        grid(false)
        plot_sif_757_c = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757_c[1:(vza_max + 1), :], levels = collect(1.2:0.2:2.6), extend = "both", cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("SIF₇₅₇ Ω = 0.4", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # NIRv Ω = 0.4
        subplot(3,3,5, polar = true)
        grid(false)
        plot_nirv_c = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), nirv_c[1:(vza_max + 1), :], levels = collect(0.20:0.03:0.44), extend = "both", cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("NIRv Ω = 0.4", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
                                            # Relative SIF 757 Ω = 0.4
        subplot(3,3,6, polar = true)
        grid(false)
        plot_sif_757_rel_c = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_757_rel_c[1:(vza_max + 1), :], levels = collect(0.0172:0.0004:0.0208), extend = "both", cmap = :viridis)
        xticks(fontsize = 11)
        yticks(fontsize = 11)
        title("Relative SIF₇₅₇ Ω = 0.4", fontsize = 14)
        colorbar(pad = 0.1).ax.tick_params(labelsize = 11)

        close(p)

    end

        #                                         # SIF Ratio
        # subplot(4,3,2, polar = true)
        # grid(false)
        # contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), sif_ratio[1:(vza_max + 1), :], cmap = :viridis)
        # xticks(fontsize = 11)
        # yticks(fontsize = 11)
        # title("SIF₇₅₇ / SIF₇₇₁")
        # colorbar(pad = 0.1).ax.tick_params(labelsize = 11)
        #                                         # LSWI
        # subplot(4,3,8, polar = true)
        # grid(false)
        # plot_ref_757 = contourf(deg2rad.(collect((0:360))), collect(0:1:vza_max), lswi[1:(vza_max + 1), :], cmap = :viridis)
        # xticks(fontsize = 11)
        # yticks(fontsize = 11)
        # title("LSWI", fontsize = 14)
        # colorbar(pad = 0.1).ax.tick_params(labelsize = 11)

    return p
end

function plot_polar_anim(dir_temp_anim::String, file_name_anim::String)

                                #### Animated GIF ####
    # Create figure for each SZA and save to a temporary directory
    if isdir(dir_temp_anim) == false # Create directory
        mkdir(dir_temp_anim)
    else
        rm(dir_temp_anim, recursive = true) # if exists, delete and recreate to avoid duplicating images
        mkdir(dir_temp_anim)
    end
    n = 1
    for i ∈ (89 - (2 * 89)):89 # SZAs
        p = plot_polar(i, vza_max, "757_vis_animated", 1.0)
        p.savefig(dir_temp_anim*"temp"*string(lpad(n, 3, '0'))*".png")
        n = n + 1
    end
    # Create an image stack for animation
    n_img = length(readdir(dir_temp_anim; join = true))
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
end