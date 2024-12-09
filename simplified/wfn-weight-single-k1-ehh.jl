# Show the contribution of the ARPES signature from a single k_1. 

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("ehh.jl")
include("wfn.jl")
include("arpes.jl")
using CairoMakie

w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])
trion = Intervalley2DChandraTrion(
    m_h = 0.21,
    m_e = 0.37,
    w = w,
    E_g = 2.84,
    E_B = 0.1,
    a = 10.3,
    b = 25.2
)

P = 1.2w
k = SA[0.1, 0.0]

kx_list = LinRange(-0.35, 0.35, 250)
k1_list = [SA[kx, 0.0] for kx in kx_list]
ω_list = LinRange(-8, 5, 200)

broaden = gaussian_broadening(10.0)
Ak1k2 = wfn(trion)

intensity_from_k1 = map(Iterators.product(k1_list, ω_list)) do (k_1, ω)
    momentum_set = momentum_calc_ehh(trion, P, k, k_1)
    k_h1 = momentum_set.k_h1
    k_h2 = momentum_set.k_h2 
    k_2 = momentum_set.k_2
    broaden(ω - E_trion(trion, P) + E_residue_ehh(trion, k_h1, k_h2)) * Ak1k2(k_1, k_2)
end

intensity_from_k1_sum = sum(intensity_from_k1, dims=1)[1, :]

f = Figure()

ax = Axis(f[1, 1])
heatmap!(ax, kx_list, ω_list, intensity_from_k1, colormap=arpes_colormap(transparency_gradience))
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

ax = Axis(f[1, 2])
lines!(ax, intensity_from_k1_sum, ω_list)
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

f
