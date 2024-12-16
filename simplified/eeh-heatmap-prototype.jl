# Plot the ARPES heatmap;
# the code hasn't been encapsulated yet.

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("eeh.jl")
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
    E_B = 0.76,
    a = 10.3,
    b = 25.2
)

P_ratio = 0.8
P = P_ratio * w

M = trion.m_h + 2trion.m_e
kx_list = LinRange(-0.35, 1.7, 250) 
k1_list = [SA[kx, 0.0] for kx in kx_list]
k1_grid = reshape(map(Iterators.product(kx_list, kx_list)) do (k_x, k_y)
    SA[k_x, k_y]
end, length(k1_list)^2)
ω_list = LinRange(-1, 3, 500) #LinRange(-8, 5, 200)

#endregion
##########################################

##########################################
#region Calculations

broaden = gaussian_broadening(20fs)
Ak1k2 = wfn(trion)

Akω_total = trion_ARPES_eeh(trion, P, k1_grid, k1_list, ω_list, Ak1k2, broaden)

#endregion
##########################################

##########################################
#region Plotting

f = Figure()
ax = Axis(f[1, 1])
heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))

# The positions of the vertical bars,
# calculated by setting k_1 = 0 or k_2 = 0.
let k_e1 = - trion.m_e / M * w + trion.m_e / M * P
    k_e2 = (trion.m_e + trion.m_h) / M * w + trion.m_e / M * P
    vlines!(ax, k_e1[1], color=colorant"lightskyblue")
    vlines!(ax, k_e2[1], color=colorant"lightskyblue")
end

# When e1 is driven out, k_1 is fixed,
# and k_2 changes freely, so we set k_2 = 0 to maximize the ARPES intensity.
dispersion_k2_zero = map(k1_list) do k
    k_2 = SA[0.0, 0.0]
    momentum_set = momentum_calc_eeh_e1(trion, P, k, k_2)
    k_h = momentum_set.k_h
    k_e2 = momentum_set.k_e2
    E_trion_eeh(trion, P) - E_residue_ehh_e1(trion, k_e2, k_h)
end
lines!(ax, kx_list, dispersion_k2_zero, color=colorant"lightskyblue")

# When e2 is driven out, k_2 is fixed,
# and k_1 changes freely, so we set k_1 = 0 to maximize the ARPES intensity.
dispersion_k1_zero = map(k1_list) do k
    k_1 = SA[0.0, 0.0]
    momentum_set = momentum_calc_eeh_e2(trion, P, k, k_1)
    k_h = momentum_set.k_h
    k_e1 = momentum_set.k_e1
    E_trion_eeh(trion, P) - E_residue_ehh_e2(trion, k_e1, k_h)
end
lines!(ax, kx_list, dispersion_k1_zero, color=colorant"lightskyblue")

ylims!(ax, (minimum(ω_list), maximum(ω_list)))
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

save("eeh-P-$(P_ratio).png", f)
f

#endregion
##########################################