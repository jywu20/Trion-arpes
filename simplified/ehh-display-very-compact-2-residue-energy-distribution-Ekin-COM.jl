# Plot the ARPES heatmap;

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("exciton.jl")
include("ehh.jl")
include("wfn.jl")
include("arpes.jl")
using CairoMakie
using LaTeXStrings
using Colors
using Distributions

##########################################
#region Parameters

caption_padding = 35

broaden = gaussian_broadening(50fs)

m_h = 0.21
m_e = 0.37
E_g = 2.67
w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])
trion = Intervalley2DChandraTrion(
    m_h = m_h,
    m_e = m_e,
    w = w,
    E_g = E_g,
    E_B = 0.75, # Binding energy for the ehh trion mode
    a = 10.3,
    b = 25.2
)
exciton = IntraValley2DExciton(
    m_h = m_h,
    m_e = m_e,
    E_g = E_g,
    E_B = 0.71,
    a = 10.4
)


# The trion momentum is set to be w
P_ratio = 1.0
P = P_ratio * w

# The exciton momentum is set to zero
Q = SA[0.0, 0.0]

# When P = w, we don't need to worry about the grid not being large enough.
# What we want is to make sure that k_1's such that k_1 = 0 or k_2 = 0
# are included in the k-grid,
# and when P=w, k_2 = 0, if and only if k_1 = -k_e.
# So for each k_e in the k-path scanned, -k_1 must be in the k-grid.
# But in this program, the k-path is a part of the k_grid,
# and it has mirror symmetry, so the condition is satisfied by default.
# When we are tilting the momentum of the trion away from P=w things can be slightly different.
kx_list = LinRange(-0.5, 0.5, 200)
ikx_Γ = argmin(abs.(kx_list))
ikx_tilted = argmin(abs.(kx_list .- -0.2))
k1_list = [SA[kx, 0.0] for kx in kx_list]
k1_grid = reshape(map(Iterators.product(kx_list, kx_list)) do (k_x, k_y)
    SA[k_x, k_y]
end, length(k1_list)^2)
ω_list = LinRange(-0.3, 3.0, 200) #LinRange(-8, 5, 200)
ω_list_plot = LinRange(-0.3, 3.0, 200)

E_c1_curve = map(k1_list) do k_h
    -E_v(exciton, k_h)
end
E_v1_curve = map(k1_list) do k_e
    E_c(exciton, k_e)
end

electron_color = :black #colorant"deepskyblue2"
hole_color = :black #colorant"coral2"
shifted_valence_color = colorant"crimson"
shifted_valence_doubled_color = colorant"darkorange"
center_line_color = colorant"gray45"
tilted_line_color = colorant"gray64"
linecut_width = 0.15
panel_width = 0.7
binding_starting_bar_pos = -0.45
binding_starting_bar_width = 0.05
binding_annotation_displacement = 0.04

σ = 0.411
μ = -0.914
σ = 0.7122748076054343
μ = -1.8521341528622541
energy_distribution = LogNormal(μ, σ)

set_theme!(fontsize=20)
f = Figure(size=(400, 300))

#endregion
##########################################

let 
    
    ##########################################
    #region Plotting the distribution of E_kin in COM frame

    ω_list = LinRange(1.5, 2.1, 200)
    # No need for the COM kinetic energy term, as Q=w => Q-w = 0.
    Ekin_COM_list = E_trion_ehh(trion, w) .- ω_list

    ax = Axis(f[1, 1], xlabel="E_kin in COM frame (eV)")
    lines!(ax, Ekin_COM_list, trion_ARPES_ehh(
        trion, w, 
        k1_grid, [SA[0.0, 0.0]], ω_list, wfn(trion),
        broaden, 0.0
    )[1, :], label=rich("Without ", rich("Z", rich("(", font=:regular), "E", rich(")", font=:regular), font=:italic), rich=:regular))
    
    lines!(ax, Ekin_COM_list, trion_ARPES_ehh(
        trion, w, 
        k1_grid, [SA[0.0, 0.0]], ω_list, wfn(trion),
        x -> pdf(energy_distribution, x),
        broaden, 0.0
    )[1, :], label=rich("With ", rich("Z", rich("(", font=:regular), "E", rich(")", font=:regular), font=:italic), rich=:regular))

    lines!(ax, Ekin_COM_list, map(x -> pdf(energy_distribution, x), Ekin_COM_list), 
        label=rich("Z", rich("(", font=:regular), "E", rich(")", font=:regular), font=:italic),
    )
    
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    axislegend(ax)
    
    save("ehh-display-very-compact-2-residue-energy-distribution-Ekin-COM.png", f)

    #endregion
    ##########################################
end
