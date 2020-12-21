using Plots

# Curves
function sensitivity_curves(model_variable::String, df::DataFrame, mat_SIF::Array, mat_REF::Array, wlf::Array, wl::Array)
    l = @layout([a{0.01h}; b c d; e f])
    title_plot = plot(title = "CliMA Top of Canopy Outgoing SIF and SW Radiation", grid = false, axis = nothing, showaxis = false, bottom_margin = -20Plots.px)

    if model_variable == "rad_ratio"
        df[!, 1] = string.(df[:, 1])
        df[!, 1] = replace.(df[:, 1], "//" => ":")
        xlab     = "Ratio Direct/Diffuse PAR"
        var_lab  = df[:, 1]
    else
        xlab    = model_variable
        var_lab = string.(model_variable, " = ", df[:, 1])
    end

    # Curves
    curve_sif757    = plot(df[:, 1], df.SIF757, ylabel          = "SIF₇₅₇ (mW m⁻² nm⁻¹ sr⁻¹)",
                           xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
    curve_sif757rel = plot(df[:, 1], df.SIF757_Relative, ylabel = "Relative SIF₇₅₇ (mW m⁻² nm⁻¹ sr⁻¹)",
                           xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
    curve_ref757    = plot(df[:, 1], df.REF757, ylabel          = "Radiance₇₅₇ (mW m⁻² nm⁻¹ sr⁻¹)",
                           xlabel = xlab, linewidth = 2, legend = false, marker = 5, framestyle = :box)
    
    # Spectra
    spectra_sif = plot()
    spectra_ref = plot()
    for i in 1:size(mat_SIF, 1)
        if i == 1
            spectra_sif = plot(wlf, mat_SIF[i, :], ylabel = "SIF₇₅₇ (mW m⁻² nm⁻¹ sr⁻¹)", xlabel = "Wavelength", linewidth = 2,
                               legend = :topright, label = var_lab[i], framestyle = :box)
        else
            spectra_sif = plot!(wlf, mat_SIF[i, :], linewidth = 2, label = var_lab[i])
        end
    end
    for i in 1:size(mat_REF, 1)
        if i == 1
            spectra_ref = plot(wl, mat_REF[i, :], ylabel = "SW Radiance (mW m⁻² nm⁻¹ sr⁻¹)", xlabel = "Wavelength", linewidth = 2,
                               legend = :topright, label = var_lab[i], framestyle = :box)
        else
            plot!(wl, mat_REF[i, :], linewidth = 2, label = var_lab[i])
        end
    end
    p = plot(title_plot, curve_sif757, curve_sif757rel, curve_ref757, spectra_sif, spectra_ref,
             plot_title = "CliMA TOC Outgoing SIF and SW Radiation", layout = l, size = (900,600))
    return p
end