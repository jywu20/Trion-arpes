# Plot the ARPES heatmap;

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("ehh.jl")
include("wfn.jl")
include("arpes.jl")
using CairoMakie
using LaTeXStrings
using Colors

##########################################
#region Parameters

w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])
trion = Intervalley2DChandraTrion(
    m_h = 0.21,
    m_e = 0.37,
    w = w,
    E_g = 2.84,
    E_B = 0.75,
    a = 10.3,
    b = 25.2
)

P_ratio = 1.0
P = P_ratio * w

M = 2trion.m_h + trion.m_e
kx_list = LinRange(-0.7, 1.7, 150) .+ trion.m_h / M * w[1] .- trion.m_h / M * P[1]
k1_list = [SA[kx, 0.0] for kx in kx_list]
k1_grid = reshape(map(Iterators.product(kx_list, kx_list)) do (k_x, k_y)
    SA[k_x, k_y]
end, length(k1_list)^2)
ω_list = LinRange(0, 3.3, 500) #LinRange(-8, 5, 200)

#endregion
##########################################

##########################################
#region Calculations

broaden = gaussian_broadening(20fs)
Ak1k2 = wfn(trion)

Akω_total = trion_ARPES_ehh(trion, P, k1_grid, k1_list, ω_list, Ak1k2, broaden)
dispersion_k_zero = map(k1_list) do k
    k_1 = SA[0.0, 0.0]
    momentum_set = momentum_calc_ehh(trion, P, k, k_1)
    k_h1 = momentum_set.k_h1
    k_h2 = momentum_set.k_h2 
    k_2  = momentum_set.k_2
    E_trion_ehh(trion, P) - E_residue_ehh(trion, k_h1, k_h2)
end

dispersion_k_equal = map(k1_list) do k
    k_1 = ((P - w) * trion.m_e / M - k) / 2
    k_2 = k_1
    momentum_set = momentum_calc_ehh(trion, P, k, k_1)
    k_h1 = momentum_set.k_h1
    k_h2 = momentum_set.k_h2 
    E_trion_ehh(trion, P) - E_residue_ehh(trion, k_h1, k_h2)
end

#endregion
##########################################

##########################################
#region Plotting

f = Figure()
ax = Axis(f[1, 1])
heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
lines!(ax, kx_list, dispersion_k_zero, label=L"k_1=0")
lines!(ax, kx_list, dispersion_k_equal, label=L"k_1=k_2")

ylims!(ax, (minimum(ω_list), maximum(ω_list)))
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
axislegend(ax)

f

#endregion
##########################################