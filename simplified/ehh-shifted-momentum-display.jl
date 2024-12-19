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

broaden = gaussian_broadening(20fs)

m_h = 0.21
m_e = 0.37
E_g = 2.84
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

# The trion momentum is set to be w
P_ratio = 1.2
P = P_ratio * w

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
ikx_tilted = argmin(abs.(kx_list .- -0.25))
k1_list = [SA[kx, 0.0] for kx in kx_list]
k1_grid = reshape(map(Iterators.product(kx_list, kx_list)) do (k_x, k_y)
    SA[k_x, k_y]
end, length(k1_list)^2)
ω_list = LinRange(0, 3.0, 200) #LinRange(-8, 5, 200)
ω_list_plot = LinRange(-0.3, 3.0, 200)

E_c1_curve = map(k1_list) do k_h
    -E_v1(trion, k_h)
end
E_v1_curve = map(k1_list) do k_e
    E_c1(trion, k_e)
end

electron_color = colorant"deepskyblue2"
hole_color = colorant"coral2"
shifted_valence_color = colorant"crimson"
shifted_valence_doubled_color = colorant"darkorange"
set_theme!(fontsize=20)
f = Figure(size=(400, 400))

#endregion
##########################################

let 
    ##########################################
    #region Calculation of trion

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
        M = 2trion.m_h + trion.m_e
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
    #region Plotting of trion

    ax = Axis(f[1, 1], xlabel="Momentum (Å⁻¹)")

    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
    lines!(ax, kx_list, E_c1_curve, color=hole_color)
    lines!(ax, kx_list, E_v1_curve, color=electron_color)
    lines!(ax, kx_list, dispersion_k_zero, label=L"k_1=0", color=shifted_valence_color, linestyle=:dash)
    lines!(ax, kx_list, dispersion_k_equal, label=L"k_1=k_2", color=shifted_valence_doubled_color, linestyle=:dash)
    
    M = m_e + 2m_h
    vlines!(ax, m_e / M * (P - w)[1], linestyle=:dot, color=:gray)
    hlines!(ax, inv_eV * m_e / 2M^2 * norm(P - w)^2 + trion.E_g - trion.E_B, linestyle=:dot, color=:gray)

    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    axislegend(ax, position=:lt)

    #endregion
    ##########################################
end

save("ehh-shifted-momentum-display.png", f)
f