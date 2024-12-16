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

kx_list = LinRange(-0.7, 1.7, 50)
ikx_Γ = argmin(abs.(kx_list))
k1_list = [SA[kx, 0.0] for kx in kx_list]
k1_grid = reshape(map(Iterators.product(kx_list, kx_list)) do (k_x, k_y)
    SA[k_x, k_y]
end, length(k1_list)^2)
ω_list = LinRange(0, 3.3, 50) #LinRange(-8, 5, 200)
ω_list_plot = LinRange(-0.3, 3.3, 50)

E_c1_curve = map(k1_list) do k_h
    -E_v(exciton, k_h)
end
E_v1_curve = map(k1_list) do k_e
    E_c(exciton, k_e)
end

electron_color = colorant"deepskyblue2"
hole_color = colorant"coral2"
shifted_valence_color = colorant"salmon1"
shifted_valence_doubled_color = colorant"gold1"
set_theme!(fontsize=20)
f = Figure(size=(1000, 400))

#endregion
##########################################

let 
    ##########################################
    #region Calculation of exciton
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
    #region Plotting the exciton
    
    ax = Axis(f[1, 1], ylabel="Energy (eV)")
    Label(f[1, 1, TopLeft()], "(a)")

    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
    lines!(ax, kx_list, dispersion, color=shifted_valence_color, linestyle=:dash)
    lines!(ax, kx_list, E_c1_curve, color=hole_color)
    lines!(ax, kx_list, E_v1_curve, color=electron_color)

    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

    #endregion
    ##########################################

    ##########################################
    #region Plotting the exciton ARPES intensity with a fixed momentum

    ax = Axis(f[1, 2], aspect=0.25, xlabel="Intensity")
    lines!(ax, Akω_total[ikx_Γ, :], ω_list)
    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax, label = false)

    #endregion
    ##########################################
end

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

    ax = Axis(f[1, 3])
    Label(f[1, 3, TopLeft()], "(b)")

    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
    lines!(ax, kx_list, E_c1_curve, color=hole_color)
    lines!(ax, kx_list, E_v1_curve, color=electron_color)
    lines!(ax, kx_list, dispersion_k_zero, label=L"k_1=0", color=shifted_valence_color, linestyle=:dash)
    lines!(ax, kx_list, dispersion_k_equal, label=L"k_1=k_2", color=shifted_valence_doubled_color, linestyle=:dash)

    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    axislegend(ax)

    #endregion
    ##########################################

    ##########################################
    #region Plotting the trion ARPES intensity with a fixed momentum

    ax = Axis(f[1, 4], aspect=0.25, xlabel="Intensity")
    lines!(ax, Akω_total[ikx_Γ, :], ω_list)
    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax, label = false)

    #endregion
    ##########################################
end

f