include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("exciton.jl")
include("wfn.jl")
include("arpes.jl")
using CairoMakie
using LaTeXStrings
using Colors

##########################################
#region Parameters

exciton = IntraValley2DExciton(
    m_h = 0.21,
    m_e = 0.37,
    E_g = 2.84,
    E_B = 0.76,
    a = 10.4,
)

Q = SA[0.0, 0.0]
M = exciton.m_e + exciton.m_h

kx_list = LinRange(-0.3, 0.3, 250)
k1_list = [SA[kx, 0.0] for kx in kx_list]
k1_grid = reshape(map(Iterators.product(kx_list, kx_list)) do (k_x, k_y)
    SA[k_x, k_y]
end, length(k1_list)^2)
ω_list = LinRange(0, 3, 500)

#endregion
##########################################

##########################################
#region Calculations

broaden = gaussian_broadening(20fs)
Ak = wfn(exciton)

Akω_total = exciton_ARPES_eh(exciton, Q, k1_list, ω_list, Ak, broaden)

#endregion
##########################################

##########################################
#region Plotting

f = Figure()
ax = Axis(f[1, 1])
heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))

ylims!(ax, (-1, 3))
xlims!(ax, (-0.5, 0.5))
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

f

#endregion
##########################################