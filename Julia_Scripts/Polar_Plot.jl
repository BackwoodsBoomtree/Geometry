include("derive_spectrum_polar.jl")

sza = 30

SIF_FR, SIF_R, reflVIS, reflNIR, wl, wlf = derive_spectrum_polar(sza)

A = reshape(reflNIR, (86,361));
B = reshape(reflVIS, (86,361));
SIFFER = reshape(SIF_R, (86,361));
SIFFER_FR = reshape(SIF_FR, (86,361));

pyplot()
PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(false)
hm = PyPlot.contourf(deg2rad.(collect((0:360))),collect(0:1:85),  A,  cmap=:viridis)
PyPlot.title("NIR reflectance BRDF")
PyPlot.colorbar()
PyPlot.gcf()

PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(false)
hm = PyPlot.contourf(deg2rad.(collect((0:360))),collect(0:1:85),  B,  cmap=:viridis)
PyPlot.title("VIS reflectance BRDF")
PyPlot.colorbar()
PyPlot.gcf()

PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(false)
hm = PyPlot.contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER, cmap=:viridis)
PyPlot.title("Red SIF emission BRDF")
PyPlot.colorbar()
PyPlot.gcf()

PyPlot.figure(figsize = (10,5))
PyPlot.subplot(1,1,1, polar=true)
PyPlot.grid(false)
hm = PyPlot.contourf(deg2rad.(collect((0:360))),collect(0:1:85),  SIFFER_FR, cmap=:viridis)
PyPlot.title("Far Red SIF emission BRDF")
PyPlot.colorbar()
PyPlot.gcf()