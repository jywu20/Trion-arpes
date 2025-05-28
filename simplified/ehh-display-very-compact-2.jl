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

caption_padding = 35

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
tilted_line_color = colorant"lightgray"
linecut_width = 0.15
panel_width = 0.7
binding_starting_bar_pos = -0.45
binding_starting_bar_width = 0.05
binding_annotation_displacement = 0.04

set_theme!(fontsize=20)
f = Figure(size=(800, 800))

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
    
    gl = f[1, 1] = GridLayout()

    ax = Axis(gl[1, 1], 
        ylabel="Energy (eV)", xticks=[-0.4, 0.0, 0.4],
        yticks=([0, exciton.E_g], ["VBM", "CBM"]),
        xtickalign = 1.0,
        ytickalign = 1.0,
    )
    Label(gl[1, 1, TopLeft()], "(a)", padding = (0, 15, 15, 0))
    Label(gl[1, 1, Bottom()], "Momentum (Å⁻¹)", padding=(0, 0, 0, caption_padding))

    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
    lines!(ax, kx_list, dispersion, color=shifted_valence_color, linestyle=:dash)
    lines!(ax, kx_list, E_c1_curve, color=hole_color)
    lines!(ax, kx_list, E_v1_curve, color=electron_color)
    vlines!(ax, kx_list[ikx_Γ], linestyle=:dot, color=center_line_color)
    hlines!(ax, exciton.E_g - exciton.E_B, linestyle=:dot, color=center_line_color)
    
    #lines!(ax, [binding_starting_bar_pos - binding_starting_bar_width, binding_starting_bar_pos + binding_starting_bar_width], [exciton.E_g, exciton.E_g], color=:black)
    arrows!(ax, [binding_starting_bar_pos], [exciton.E_g], [0.0], [-0.9exciton.E_B], color=:black)
    text!(ax, binding_starting_bar_pos + binding_annotation_displacement, exciton.E_g - exciton.E_B/2, 
        text=rich("E", subscript("B,X", font=:regular), font=:italic)
    ) 

    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

    #endregion
    ##########################################

    ##########################################
    #region Plotting the exciton ARPES intensity with a fixed momentum

    ax = Axis(gl[1, 2], xticks=([0.0], [" "]), xtickcolor=:white)
    # Labels of line cuts are removed.
    #Label(gl[1, 2, Top()], "(b)", padding = (0, 0, 15, 0), halign=:left)
    Label(gl[1, 2, Bottom()], "Intensity", padding=(0, 0, 0, caption_padding))

    # Usually we can just hide x decorations,
    # but here to align the label "Intensity" with the label "Momentum",
    # we have to make the x ticks and the x tick labels invisible
    # so that they still occupy their space.
    #hidexdecorations!(ax, label = false)
    
    lines!(ax, Akω_total[ikx_Γ, :], ω_list, color=center_line_color)
    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidexdecorations!(ax, ticks = false, ticklabels = false, label = false)
    hideydecorations!(ax)
    # Note that in a Makie plot there are two aspect ratios:
    # the first is the aspect ratio of a cell in the Figure,
    # and the second is the aspect ratio of an Axis in the cell.
    # If we set the latter, then unwanted empty spaces will appear around the Axis.
    # Here we set the width of column 2 to be 0.25 times the height of row 1,
    # which means we are setting the aspect ratio of the cell,
    # and it will automatically determine the aspect ratio of the Axis within the cell.
    colsize!(gl, 1, Aspect(1, panel_width))
    colsize!(gl, 2, Aspect(1, linecut_width))
    #colsize!(f.layout, 1, Aspect(1, 1.2))
    colgap!(gl, 0)

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

    gl = f[1, 2] = GridLayout()

    ax = Axis(gl[1, 1], xticks=[-0.4, 0.0, 0.4], 
        yticks=([0, exciton.E_g], ["", ""]),
        xtickalign = 1.0,
        ytickalign = 1.0,
    )
    colsize!(gl, 1, Aspect(1, panel_width))
    # Because labels of line cuts are removed, relabeling is needed.
    Label(gl[1, 1, TopLeft()], "(b)", padding = (0, 15, 15, 0))
    Label(gl[1, 1, Bottom()], "Momentum (Å⁻¹)", padding=(0, 0, 0, caption_padding))

    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
    lines!(ax, kx_list, E_c1_curve, color=hole_color)
    lines!(ax, kx_list, E_v1_curve, color=electron_color)
    lines!(ax, kx_list, dispersion_k_zero, label=L"k_1=0", color=shifted_valence_color, linestyle=:dash)
    lines!(ax, kx_list, dispersion_k_equal, label=L"k_1=k_2", color=shifted_valence_doubled_color, linestyle=:dash)
    vlines!(ax, kx_list[ikx_Γ], linestyle=:dot, color=center_line_color)
    vlines!(ax, kx_list[ikx_tilted], linestyle=:dot, color=tilted_line_color)
    hlines!(ax, trion.E_g - trion.E_B, linestyle=:dot, color=center_line_color)

    arrows!(ax, [binding_starting_bar_pos], [trion.E_g], [0.0], [-0.9trion.E_B], color=:black)
    text!(ax, binding_starting_bar_pos + binding_annotation_displacement, trion.E_g - trion.E_B/2, text=rich("E", subscript("B,T", font=:regular), font=:italic)) 

    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    axislegend(ax)

    #endregion
    ##########################################

    ##########################################
    #region Plotting the trion ARPES intensity with a fixed momentum

    ax = Axis(gl[1, 2], xticks=([0.0], [" "]), xtickcolor=:white,)
    hideydecorations!(ax)
    # Labels of line cuts are removed.
    #Label(gl[1, 3, Top()], "(e)", padding = (0, 15, 15, 0), halign=:left)
    
    lines!(ax, Akω_total[ikx_tilted, :], ω_list, color=tilted_line_color)
    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    #hidexdecorations!(ax, label = false)
    colsize!(gl, 2, Aspect(1, linecut_width))

    # The reason not to use hidexdecorations is discussed above 
    ax = Axis(gl[1, 3], xticks=([0.0], [" "]), xtickcolor=:white)
    hideydecorations!(ax)
    # Labels of line cuts are removed.
    #Label(gl[1, 2, Top()], "(d)", padding = (0, 0, 15, 0), halign=:left)

    lines!(ax, Akω_total[ikx_Γ, :], ω_list, color=center_line_color)
    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    #hidexdecorations!(ax, label = false)
    colsize!(gl, 3, Aspect(1, linecut_width))

    colgap!(gl, 0)

    Label(gl[1, 2:3, Bottom()], "Intensity", tellwidth=false, tellheight=false, padding=(0, 0, 0, caption_padding))

    #endregion
    ##########################################
end

P_ratio = 1.2
Q = (P_ratio - 1) * w

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
    
    gl = f[2, 1] = GridLayout()

    ax = Axis(gl[1, 1], ylabel="Energy (eV)", xticks=[-0.4, 0.0, 0.4],
        yticks=([0, exciton.E_g], ["VBM", "CBM"]),
        xtickalign = 1.0,
        ytickalign = 1.0,
    )
    Label(gl[1, 1, TopLeft()], "(c)", padding = (0, 15, 15, 0))
    Label(gl[1, 1, Bottom()], "Momentum (Å⁻¹)", padding=(0, 0, 0, caption_padding))

    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
    lines!(ax, kx_list, dispersion, color=shifted_valence_color, linestyle=:dash)
    lines!(ax, kx_list, E_c1_curve, color=hole_color)
    lines!(ax, kx_list, E_v1_curve, color=electron_color)
    
    M = m_e + m_h
    vlines!(ax, m_e / M * Q[1], linestyle=:dot, color=center_line_color)
    hlines!(ax, inv_eV * m_e / 2M^2 * norm(Q)^2 + exciton.E_g - exciton.E_B, linestyle=:dot, color=center_line_color)
    
    #arrows!(ax, [binding_starting_bar_pos], [exciton.E_g], [0.0], [-0.9exciton.E_B], color=:black)
    #text!(ax, binding_starting_bar_pos + binding_annotation_displacement, exciton.E_g - exciton.E_B/2, 
    #    text=rich("E", subscript("B,X", font=:regular), font=:italic)
    #) 

    ikx_center = argmin(abs.(kx_list .- m_e / M * Q[1]))

    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

    #endregion
    ##########################################

    ##########################################
    #region Plotting the exciton ARPES intensity with a fixed momentum

    ax = Axis(gl[1, 2], xticks=([0.0], [" "]), )
    # Labels of line cuts are removed.
    #Label(gl[1, 2, Top()], "(b)", padding = (0, 0, 15, 0), halign=:left)
    #Label(gl[1, 2, Bottom()], "Intensity", padding=(0, 0, 0, caption_padding))

    # Usually we can just hide x decorations,
    # but here to align the label "Intensity" with the label "Momentum",
    # we have to make the x ticks and the x tick labels invisible
    # so that they still occupy their space.
    #hidexdecorations!(ax, label = false)
    
    lines!(ax, Akω_total[ikx_center, :], ω_list, color=center_line_color)
    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidexdecorations!(ax, ticks = false, ticklabels = false, label = false)
    hideydecorations!(ax)
    # Note that in a Makie plot there are two aspect ratios:
    # the first is the aspect ratio of a cell in the Figure,
    # and the second is the aspect ratio of an Axis in the cell.
    # If we set the latter, then unwanted empty spaces will appear around the Axis.
    # Here we set the width of column 2 to be 0.25 times the height of row 1,
    # which means we are setting the aspect ratio of the cell,
    # and it will automatically determine the aspect ratio of the Axis within the cell.
    colsize!(gl, 1, Aspect(1, panel_width))
    colsize!(gl, 2, Aspect(1, linecut_width))
    #colsize!(f.layout, 1, Aspect(1, 1.2))
    colgap!(gl, 0)

    #endregion
    ##########################################
end

# The trion momentum is set to be w plus the displacement
P = P_ratio * w

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

    gl = f[2, 2] = GridLayout()
    ax = Axis(gl[1, 1], xticks=[-0.4, 0.0, 0.4], 
        yticks=([0, exciton.E_g], ["", ""]),
        xtickalign = 1.0,
        ytickalign = 1.0,
    )
    colsize!(gl, 1, Aspect(1, panel_width))

    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
    lines!(ax, kx_list, E_c1_curve, color=hole_color)
    lines!(ax, kx_list, E_v1_curve, color=electron_color)
    lines!(ax, kx_list, dispersion_k_zero, label=L"k_1=0", color=shifted_valence_color, linestyle=:dash)
    lines!(ax, kx_list, dispersion_k_equal, label=L"k_1=k_2", color=shifted_valence_doubled_color, linestyle=:dash)
    
    M = m_e + 2m_h
    hlines!(ax, inv_eV * m_e / 2M^2 * norm(P - w)^2 + trion.E_g - trion.E_B, linestyle=:dot, color=center_line_color)

    ik_center = argmin(abs.(kx_list .- m_e / M * (P - w)[1]))
    ikx_tilted = argmin(abs.(kx_list .- (m_e / M * (P - w)[1] -0.25)))
    vlines!(ax, m_e / M * (P - w)[1], linestyle=:dot, color=center_line_color)
    vlines!(ax, kx_list[ikx_tilted], linestyle=:dot, color=tilted_line_color)
    
    #arrows!(ax, [binding_starting_bar_pos], [trion.E_g], [0.0], [-0.9trion.E_B], color=:black)
    #text!(ax, binding_starting_bar_pos + binding_annotation_displacement, trion.E_g - trion.E_B/2, text=rich("E", subscript("B,T", font=:regular), font=:italic)) 

    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    axislegend(ax, position=:rb)
    #colsize!(f.layout, 3, Aspect(1, panel_width))
    
    # Because labels of line cuts are removed, relabeling is needed.
    Label(gl[1, 1, TopLeft()], "(d)", padding = (0, 15, 15, 0))
    Label(gl[1, 1, Bottom()], "Momentum (Å⁻¹)", padding=(0, 0, 0, caption_padding))

    # The statements below are for keeping panel (d) in the correct position.

    ax = Axis(gl[1, 2], xticks=([0.0], [" "]), xtickcolor=:white,)
    lines!(ax, Akω_total[ikx_tilted, :], ω_list, color=tilted_line_color)
    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hideydecorations!(ax)
    colsize!(gl, 2, Aspect(1, linecut_width))

    ax = Axis(gl[1, 3], xticks=([0.0], [" "]), xtickcolor=:white,)
    lines!(ax, Akω_total[ik_center, :], ω_list, color=center_line_color)
    ylims!(ax, (minimum(ω_list_plot), maximum(ω_list_plot)))
    hideydecorations!(ax)
    colsize!(gl, 3, Aspect(1, linecut_width))
    colgap!(gl, 0)


    #endregion
    ##########################################
end

#rowsize!(f.layout, 1, Fixed(200))
#rowsize!(f.layout, 2, Fixed(200))

save("ehh-display-very-compact-2.png", f)
f