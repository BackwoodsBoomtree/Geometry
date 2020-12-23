using Plots

# Curves
function sensitivity_curves(model_variable::String, df::DataFrame, mat_SIF::Array, mat_REF::Array, wlf::Array, wl::Array)

    if model_variable == "rad_ratio"
        df[!, 1] = string.(df[:, 1])
        df[!, 1] = replace.(df[:, 1], "//" => ":")
        xlab     = "Ratio Direct/Diffuse PAR"
        var_lab  = df[:, 1]
    elseif model_variable == "cab_gradient"
        xlab = "Cab Gradient"
        var_lab = df[:, 1]
    elseif model_variable == "cab_gradient_year"
        xlab = "Cab Gradient Ying et al. 2017"
        var_lab = df[:, 1]
    else
        xlab    = model_variable
        var_lab = string.(model_variable, " = ", df[:, 1])
    end

    if mat_SIF == mat_SIF_sunlit || mat_SIF == mat_SIF_shaded || mat_SIF == mat_SIF_scattered || mat_SIF == mat_SIF_soil
        l = @layout([a{0.01h}; b; c d])
        if mat_SIF == mat_SIF_sunlit
            plot_title = "CliMA Top of Canopy Outgoing SIF and Reflectance (mW m⁻² nm⁻¹ sr⁻¹) - SIF Sunlit"
            curve_sif757    = plot(df[:, 1], df.SIF757_sunlit, ylabel          = "SIF₇₅₇",
                                xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        elseif mat_SIF == mat_SIF_shaded
            plot_title = "CliMA Top of Canopy Outgoing SIF and Reflectance (mW m⁻² nm⁻¹ sr⁻¹) - SIF Shaded"
            curve_sif757    = plot(df[:, 1], df.SIF757_shaded, ylabel          = "SIF₇₅₇",
                                xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        elseif mat_SIF == mat_SIF_scattered
            plot_title = "CliMA Top of Canopy Outgoing SIF and Reflectance (mW m⁻² nm⁻¹ sr⁻¹) - SIF Scattered"
            curve_sif757    = plot(df[:, 1], df.SIF757_scattered, ylabel          = "SIF₇₅₇",
                                xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        elseif mat_SIF == mat_SIF_soil
            plot_title = "CliMA Top of Canopy Outgoing SIF and Reflectance (mW m⁻² nm⁻¹ sr⁻¹) - SIF Soil"
            curve_sif757    = plot(df[:, 1], df.SIF757_soil, ylabel          = "SIF₇₅₇",
                                xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        end
    else
        l = @layout([a{0.01h}; b c d; e f])
        plot_title = "CliMA Top of Canopy Outgoing SIF and Reflectance (mW m⁻² nm⁻¹ sr⁻¹)"
        # Curves
        curve_sif757    = plot(df[:, 1], df.SIF757, ylabel          = "SIF₇₅₇",
                               xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        curve_sif757rel = plot(df[:, 1], df.SIF757_Relative, ylabel = "Relative SIF₇₅₇",
                               xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
        curve_ref757    = plot(df[:, 1], df.REF757, ylabel          = "Reflectance₇₅₇",
                               xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
    end
    
    # Spectra Plots
    spectra_sif = plot()
    spectra_ref = plot()

    # SIF Spectrum
    for i in 1:size(mat_SIF, 1)
        if i == 1
            spectra_sif = plot(wlf, mat_SIF[i, :], ylabel = "SIF₇₅₇", xlabel = "Wavelength", linewidth = 2,
                               legend = :topright, label = var_lab[i], framestyle = :box)
        else
            spectra_sif = plot!(wlf, mat_SIF[i, :], linewidth = 2, label = var_lab[i])
        end
    end

    # Reflectance Spectrum
    for i in 1:size(mat_REF, 1)
        if i == 1
            spectra_ref = plot(wl, mat_REF[i, :], ylabel = "Reflectance", xlabel = "Wavelength", linewidth = 2,
                               legend = :topright, label = var_lab[i], framestyle = :box)
        else
            plot!(wl, mat_REF[i, :], linewidth = 2, label = var_lab[i])
        end
    end
    
    # Title Plot
    title_plot = plot(title = plot_title, grid = false, axis = nothing, showaxis = false, bottom_margin = -20Plots.px)

    # Final Plot
    if mat_SIF == mat_SIF_sunlit || mat_SIF == mat_SIF_shaded || mat_SIF == mat_SIF_scattered || mat_SIF == mat_SIF_soil
        p = plot(title_plot, curve_sif757, spectra_sif, spectra_ref, layout = l, size = (900,600))
    else
        p = plot(title_plot, curve_sif757, curve_sif757rel, curve_ref757, spectra_sif, spectra_ref, layout = l, size = (900,600))
    end

    return p
end