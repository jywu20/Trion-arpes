# Show the contribution of the ARPES signature from a single k_1. 

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("ehh.jl")
include("wfn.jl")
include("arpes.jl")
using CairoMakie
using Colors

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
k = SA[0.36, 0.0]

M = 2trion.m_h + trion.m_e
kx_list = LinRange(-0.7, 1.7, 150) .+ trion.m_h / M * w[1] .- trion.m_h / M * P[1]
k1_list = [SA[kx, 0.0] for kx in kx_list]
k1_grid = reshape(map(Iterators.product(kx_list, kx_list)) do (k_x, k_y)
    SA[k_x, k_y]
end, length(k1_list)^2)
ω_list = LinRange(0, 3, 500) #LinRange(-8, 5, 200)

broaden = gaussian_broadening(20fs)
Ak1k2 = wfn(trion)

intensity_from_k1 = map(Iterators.product(k1_list, ω_list)) do (k_1, ω)
    momentum_set = momentum_calc_ehh(trion, P, k, k_1)
    k_h1 = momentum_set.k_h1
    k_h2 = momentum_set.k_h2 
    k_2 = momentum_set.k_2
    broaden(ω - E_trion_ehh(trion, P) + E_residue_ehh(trion, k_h1, k_h2)) * abs(Ak1k2(k_1, k_2))^2
end

intensity_from_k1_sum = sum(intensity_from_k1, dims=1)[1, :]

dispersion_k1 = map(k1_list) do k_1
    momentum_set = momentum_calc_ehh(trion, P, k, k_1)
    k_h1 = momentum_set.k_h1
    k_h2 = momentum_set.k_h2 
    E_trion_ehh(trion, P) - E_residue_ehh(trion, k_h1, k_h2)
end

Akω_total = trion_ARPES_ehh(trion, P, k1_grid, k1_list, ω_list, Ak1k2, broaden)

##

f = Figure(size=(1000,500))

ax = Axis(f[1, 1])
lines!(ax, intensity_from_k1_sum, ω_list)
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

ax = Axis(f[1, 2])
lines!(ax, kx_list, dispersion_k1, color=colorant"lightskyblue")
heatmap!(ax, kx_list, ω_list, intensity_from_k1, colormap=arpes_colormap(transparency_gradience))
ylims!(ax, (minimum(ω_list), maximum(ω_list)))
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

ax = Axis(f[1, 3])
heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

f
