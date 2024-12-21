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

M = trion.m_h + 2trion.m_e
kx_list = LinRange(-0.35, 1.7, 250) 
k1_list = [SA[kx, 0.0] for kx in kx_list]
k1_grid = reshape(map(Iterators.product(kx_list, kx_list)) do (k_x, k_y)
    SA[k_x, k_y]
end, length(k1_list)^2)
ω_list = LinRange(-0.3, 3.0, 200)

electron_color = colorant"deepskyblue2"
hole_color = colorant"coral2"
shifted_valence_color = colorant"crimson"
shifted_effective_mass_color = colorant"darkorange"
set_theme!(fontsize=20)
f = Figure(size=(1000, 400))

E_c1_curve = map(k1_list) do k_h
    -E_v1(trion, k_h)
end
E_v1_curve = map(k1_list) do k_e
    E_c1(trion, k_e)
end
E_c2_curve = map(k1_list) do k_h
    -E_v2(trion, k_h)
end
E_v2_curve = map(k1_list) do k_e
    E_c2(trion, k_e)
end

#endregion
##########################################

# The zero-momentum case
let 
    P_ratio = 1.0
    P = P_ratio * w
    
    ##########################################
    #region Calculations

    broaden = gaussian_broadening(20fs)
    Ak1k2 = wfn(trion)

    Akω_total = trion_ARPES_eeh(trion, P, k1_grid, k1_list, ω_list, Ak1k2, broaden)

    #endregion
    ##########################################

    ##########################################
    #region Plotting
    
    ax = Axis(f[1, 1], ylabel="Energy (eV)", xlabel="Momentum (Å⁻¹)")
    Label(f[1, 1, TopLeft()], "(a)", padding = (0, 15, 15, 0))
    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))

    ylims!(ax, (minimum(ω_list), maximum(ω_list)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

    #endregion
    ##########################################
    
    ##########################################
    #region Labeling

    # The positions of the vertical bars,
    # calculated by setting k_1 = 0 or k_2 = 0.
    let k_e1 = - trion.m_e / M * w + trion.m_e / M * P
        k_e2 = (trion.m_e + trion.m_h) / M * w + trion.m_e / M * P
        vlines!(ax, k_e1[1], linestyle=:dot, color=:gray)
        vlines!(ax, k_e2[1], linestyle=:dot, color=:gray)
    end

    # When e1 is driven out, k_1 is fixed,
    # and k_2 changes freely, so we set k_2 = 0 to maximize the ARPES intensity.
    dispersion_k2_zero = map(k1_list) do k
        k_2 = SA[0.0, 0.0]
        momentum_set = momentum_calc_eeh_e1(trion, P, k, k_2)
        k_h = momentum_set.k_h
        k_e2 = momentum_set.k_e2
        E_trion_eeh(trion, P) - E_residue_eeh_e1(trion, k_e2, k_h)
    end
    lines!(ax, kx_list, dispersion_k2_zero, color=shifted_valence_color, linestyle=:dash)

    # When e2 is driven out, k_2 is fixed,
    # and k_1 changes freely, so we set k_1 = 0 to maximize the ARPES intensity.
    dispersion_k1_zero = map(k1_list) do k
        k_1 = SA[0.0, 0.0]
        momentum_set = momentum_calc_eeh_e2(trion, P, k, k_1)
        k_h = momentum_set.k_h
        k_e1 = momentum_set.k_e1
        E_trion_eeh(trion, P) - E_residue_eeh_e2(trion, k_e1, k_h)
    end
    lines!(ax, kx_list, dispersion_k1_zero, color=shifted_valence_color, linestyle=:dash)

    # The maximal frequency that can be achieved
    dispersion_k1_max = map(k1_list) do k
        E_trion_eeh(trion, P) - trion.E_g - inv_eV * norm(P - w - k)^2 / 2(trion.m_h + trion.m_e)
    end
    dispersion_k2_max = map(k1_list) do k
        E_trion_eeh(trion, P) - trion.E_g - inv_eV * norm(P - k)^2 / 2(trion.m_h + trion.m_e)
    end
    lines!(ax, kx_list, dispersion_k1_max, color=shifted_effective_mass_color, linestyle=:dash)
    lines!(ax, kx_list, dispersion_k2_max, color=shifted_effective_mass_color, linestyle=:dash)

    lines!(ax, kx_list, E_c1_curve, color=electron_color)
    lines!(ax, kx_list, E_c2_curve, color=electron_color)
    lines!(ax, kx_list, E_v1_curve, color=hole_color)
    lines!(ax, kx_list, E_v2_curve, color=hole_color)

    #endregion
    ##########################################
end

let 
    P_ratio = 1.2
    P = P_ratio * w
    
    ##########################################
    #region Calculations

    broaden = gaussian_broadening(20fs)
    Ak1k2 = wfn(trion)

    Akω_total = trion_ARPES_eeh(trion, P, k1_grid, k1_list, ω_list, Ak1k2, broaden)

    #endregion
    ##########################################

    ##########################################
    #region Plotting
    
    ax = Axis(f[1, 2], ylabel="Energy (eV)", xlabel="Momentum (Å⁻¹)")
    Label(f[1, 2, TopLeft()], "(b)", padding = (0, 15, 15, 0))
    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))

    ylims!(ax, (minimum(ω_list), maximum(ω_list)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

    #endregion
    ##########################################
    
    ##########################################
    #region Labeling

    # The positions of the vertical bars,
    # calculated by setting k_1 = 0 or k_2 = 0.
    let k_e1 = - trion.m_e / M * w + trion.m_e / M * P
        k_e2 = (trion.m_e + trion.m_h) / M * w + trion.m_e / M * P
        vlines!(ax, k_e1[1], linestyle=:dot, color=:gray)
        vlines!(ax, k_e2[1], linestyle=:dot, color=:gray)
    end

    # When e1 is driven out, k_1 is fixed,
    # and k_2 changes freely, so we set k_2 = 0 to maximize the ARPES intensity.
    dispersion_k2_zero = map(k1_list) do k
        k_2 = SA[0.0, 0.0]
        momentum_set = momentum_calc_eeh_e1(trion, P, k, k_2)
        k_h = momentum_set.k_h
        k_e2 = momentum_set.k_e2
        E_trion_eeh(trion, P) - E_residue_eeh_e1(trion, k_e2, k_h)
    end
    lines!(ax, kx_list, dispersion_k2_zero, color=shifted_valence_color, linestyle=:dash)

    # When e2 is driven out, k_2 is fixed,
    # and k_1 changes freely, so we set k_1 = 0 to maximize the ARPES intensity.
    dispersion_k1_zero = map(k1_list) do k
        k_1 = SA[0.0, 0.0]
        momentum_set = momentum_calc_eeh_e2(trion, P, k, k_1)
        k_h = momentum_set.k_h
        k_e1 = momentum_set.k_e1
        E_trion_eeh(trion, P) - E_residue_eeh_e2(trion, k_e1, k_h)
    end
    lines!(ax, kx_list, dispersion_k1_zero, color=shifted_valence_color, linestyle=:dash)

    # The maximal frequency that can be achieved
    dispersion_k1_max = map(k1_list) do k
        E_trion_eeh(trion, P) - trion.E_g - inv_eV * norm(P - w - k)^2 / 2(trion.m_h + trion.m_e)
    end
    dispersion_k2_max = map(k1_list) do k
        E_trion_eeh(trion, P) - trion.E_g - inv_eV * norm(P - k)^2 / 2(trion.m_h + trion.m_e)
    end
    lines!(ax, kx_list, dispersion_k1_max, color=shifted_effective_mass_color, linestyle=:dash)
    lines!(ax, kx_list, dispersion_k2_max, color=shifted_effective_mass_color, linestyle=:dash)

    lines!(ax, kx_list, E_c1_curve, color=electron_color)
    lines!(ax, kx_list, E_c2_curve, color=electron_color)
    lines!(ax, kx_list, E_v1_curve, color=hole_color)
    lines!(ax, kx_list, E_v2_curve, color=hole_color)

    #endregion
    ##########################################
end

save("eeh-display.png", f)
f