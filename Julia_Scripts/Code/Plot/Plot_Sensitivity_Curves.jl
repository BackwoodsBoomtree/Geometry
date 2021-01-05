using Plots

# Curves
function sensitivity_curves(model_variable::String, df::DataFrame, mat_SIF::Array, mat_reflectance::Array, wlf::Array, wl::Array)

    plot_title = "CliMA Top of Canopy SIF, Radiance, Reflectance, and VIs (mW m⁻² nm⁻¹ sr⁻¹)"

    if model_variable == "rad_ratio"
        df[!, 1] = string.(df[:, 1])
        df[!, 1] = replace.(df[:, 1], "//" => ":")
        xlab     = "Ratio Direct/Diffuse PAR"
        var_lab  = df[:, 1]
    elseif model_variable == "cab_gradient"
        xlab = "Cab Gradient"
        var_lab = df[:, 1]
    elseif model_variable == "cab_Yang_2017"
        xlab = "Cab Gradient Yang et al. 2017"
        var_lab = df[:, 1]
    else
        xlab    = model_variable
        var_lab = string.(model_variable, " = ", df[:, 1])
    end

    if mat_SIF == mat_SIF_obs_sunlit || mat_SIF == mat_SIF_obs_shaded || mat_SIF == mat_SIF_obs_scattered || mat_SIF == mat_SIF_obs_soil
        if mat_SIF == mat_SIF_obs_sunlit
            plot_title = "CliMA Top of Canopy SIF, Radiance, Reflectance, and VIs (mW m⁻² nm⁻¹ sr⁻¹) - SIF Sunlit"
            curve_sif757    = plot(df[:, 1], df.SIF757_sunlit, ylabel        = "SIF₇₅₇",
                                   xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        elseif mat_SIF == mat_SIF_obs_shaded
            plot_title = "CliMA Top of Canopy SIF, Radiance, Reflectance, and VIs (mW m⁻² nm⁻¹ sr⁻¹) - SIF Shaded"
            curve_sif757    = plot(df[:, 1], df.SIF757_shaded, ylabel        = "SIF₇₅₇",
                                   xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        elseif mat_SIF == mat_SIF_obs_scattered
            plot_title = "CliMA Top of Canopy SIF, Radiance, Reflectance, and VIs (mW m⁻² nm⁻¹ sr⁻¹) - SIF Scattered"
            curve_sif757    = plot(df[:, 1], df.SIF757_scattered, ylabel     = "SIF₇₅₇",
                                   xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        elseif mat_SIF == mat_SIF_obs_soil
            plot_title = "CliMA Top of Canopy SIF, Radiance, Reflectance, and VIs (mW m⁻² nm⁻¹ sr⁻¹) - SIF Soil"
            curve_sif757    = plot(df[:, 1], df.SIF757_soil, ylabel          = "SIF₇₅₇",
                                   xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        end
        curve_plots             = plot(curve_sif757)
    else
        # Curves
        curve_layout            = @layout([a b c])
        curve_sif757            = plot(df[:, 1], df.SIF757, ylabel          = "SIF₇₅₇",
                                        xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        curve_sif757rel         = plot(df[:, 1], df.SIF757_Relative, ylabel = "Relative SIF₇₅₇",
                                        xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        curve_reflectance757    = plot(df[:, 1], df.Rad757, ylabel          = "Radiance₇₅₇",
                                        xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        curve_plots             = plot(curve_sif757, curve_sif757rel, curve_reflectance757, layout = curve_layout)
    end
    
    # Spectra Plots
    spectra_layout         = @layout([a b])
    spectra_sif            = plot()
    spectra_reflectance    = plot()
    for i in 1:size(mat_SIF, 1)                   # SIF Spectrum
        if i == 1
            spectra_sif = plot(wlf, mat_SIF[i, :], ylabel = "SIF", xlabel = "Wavelength", linewidth = 2,
                               legend = :topright, label = var_lab[i], framestyle = :box)
        else
            spectra_sif = plot!(wlf, mat_SIF[i, :], linewidth = 2, label = var_lab[i])
        end
    end
    for i in 1:size(mat_reflectance, 1)          # Reflectance Spectrum
        if i == 1
            spectra_reflectance = plot(wl, mat_reflectance[i, :], ylabel = "Reflectance", xlabel = "Wavelength", linewidth = 2,
                               legend = :topright, label = var_lab[i], framestyle = :box)
        else
            plot!(wl, mat_reflectance[i, :], linewidth = 2, label = var_lab[i])
        end
    end
    if mat_SIF == mat_SIF_obs_sunlit || mat_SIF == mat_SIF_obs_shaded || mat_SIF == mat_SIF_obs_scattered || mat_SIF == mat_SIF_obs_soil
        spectra_plots = plot(spectra_sif)
    else
        spectra_plots = plot(spectra_sif, spectra_reflectance, layout = spectra_layout)
    end

    # VI Plots
    vi_layout    = @layout([a b c d; e{0.01h}])
    curve_NDVI   = plot(df[:, 1], df.NDVI, ylabel         = "NDVI",
                            linewidth = 2, legend = false, marker = 5, framestyle = :box)
    curve_NIRv   = plot(df[:, 1], df.NIRv, ylabel         = "NIRv",
                            linewidth = 2, legend = false, marker = 5, framestyle = :box)
    curve_EVI    = plot(df[:, 1], df.EVI, ylabel          = "EVI",
                            linewidth = 2, legend = false, marker = 5, framestyle = :box)
    curve_LSWI   = plot(df[:, 1], df.LSWI, ylabel         = "LSWI",
                            linewidth = 2, legend = false, marker = 5, framestyle = :box)
    # VI Label
    xlab_VIs     = plot(title = xlab, grid = false, axis = nothing, showaxis = false, titlefontsize = 11, bottom_margin = -20Plots.px)
    vi_plots     = plot(curve_NDVI, curve_NIRv, curve_EVI, curve_LSWI, xlab_VIs, layout = vi_layout, bottom_margin = -10Plots.px)

    # Title Plot
    title_plot = plot(title = plot_title, grid = false, axis = nothing, showaxis = false, bottom_margin = -30Plots.px)

    # Final Plot
    if mat_SIF == mat_SIF_obs_sunlit || mat_SIF == mat_SIF_obs_shaded || mat_SIF == mat_SIF_obs_scattered || mat_SIF == mat_SIF_obs_soil
        l = @layout([a{0.01h}; b; c])
        p = plot(title_plot, curve_plots, spectra_plots, layout = l, size = (900,600))
    else
        l = @layout([a{0.01h}; b; c; d])
        p = plot(title_plot, curve_plots, spectra_plots, vi_plots, layout = l, size = (900,750))
    end

    return p
end