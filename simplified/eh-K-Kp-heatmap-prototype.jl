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

w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])
exciton = InterValley2DExciton(
    m_h = 0.21,
    m_e = 0.37,
    E_g = 2.84,
    E_B = 0.76,
    a = 10.4,
    w = w
)

Q_ratio = 0.8
Q = Q_ratio * w
M = exciton.m_e + exciton.m_h

kx_list = LinRange(-0.35, 1.7, 250) 
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
dispersion = map(k1_list) do k
    k_e = k
    momentum_set = momentum_calc(exciton, Q, k)
    k_h = momentum_set.k_h
    E_exciton(exciton, Q) - E_v(exciton, k_h)
end

#endregion
##########################################

##########################################
#region Plotting

f = Figure()
ax = Axis(f[1, 1])
heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
lines!(ax, kx_list, dispersion, color=colorant"lightblue1")

ylims!(ax, (-1, 3))
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

save("eh-K-Kp-$Q_ratio.png", f)

f

#endregion
##########################################