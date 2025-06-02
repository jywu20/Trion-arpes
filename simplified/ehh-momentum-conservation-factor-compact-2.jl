# Plot the contribution of the trion wave function structural factor and the energy conservation delta function

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("exciton.jl")
include("ehh.jl")
include("eeh.jl")
include("wfn.jl")
include("arpes.jl")
using CairoMakie
using LaTeXStrings
using Colors

function cross_mark(ax::Axis, pos::AbstractVector, len::Float64)
    lines!(ax, [pos[1] - len, pos[1] + len], [pos[2], pos[2]], color=:white, linewidth=0.5)
    lines!(ax, [pos[1], pos[1]], [pos[2] - len, pos[2] + len], color=:white, linewidth=0.5)
end

##########################################
#region Parameters

set_theme!(fontsize=20)
label_shift = 0.05
f = Figure(size=(600, 850), figure_padding=0.0)

broaden = gaussian_broadening(20fs)

m_h = 0.21
m_e = 0.37
E_g = 2.84
w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])

kx_list = LinRange(-0.3, 0.2, 500)
ky_list = LinRange(-0.2, 0.3, 500)
k1_grid = map(Iterators.product(kx_list, ky_list)) do (k_x, k_y)
    SA[k_x, k_y]
end

let 
    trion = Intervalley2DChandraTrion(
        m_h = m_h,
        m_e = m_e,
        w = w,
        E_g = E_g,
        E_B = 0.75, # Binding energy for the ehh trion mode
        a = 10.3,
        b = 25.2
    )
    M = 2m_h + m_e
    A_k1k2 = wfn(trion)
    
    k = SVector{2, Float64}([0.3, -0.1])
    P = 1.2w + SA[0.0, 0.1]
    
    #endregion
    ##########################################
    
    # The k_1=0 case
    let 
        k_1_max, k_2_max, ω = let k_1 = SA[0.0, 0.0]
            momenta = momentum_calc_ehh(trion, P, k, SA[0.0, 0.0])
            k_h1 = momenta.k_h1
            k_h2 = momenta.k_h2
            ω = E_trion_ehh(trion, P) - E_residue_ehh(trion, k_h1, k_h2)
            k_2 = momenta.k_2
            k_1, k_2, ω
        end
    
        ##########################################
        #region Calculating the k_1 = 0 case
        
        E_residue = map(k1_grid) do k_1
            k_h1 = momentum_calc_ehh(trion, P, k, k_1).k_h1
            k_h2 = momentum_calc_ehh(trion, P, k, k_1).k_h2 
            broaden(ω - E_trion_ehh(trion, P) + E_residue_ehh(trion, k_h1, k_h2))
        end
        
        Asq = map(k1_grid) do k_1
            k_2 = momentum_calc_ehh(trion, P, k, k_1).k_2
            abs2(A_k1k2(k_1, k_2))
        end
        
        #endregion
        ##########################################
        
        ##########################################
        #region Plotting the k_1 = 0 case
        
        ax = Axis(f[1, 1], title=L"|A|^2", titlefont=:regular, aspect=1)
        heatmap!(ax, kx_list, ky_list, Asq,
            colorrange=(minimum(Asq), maximum(Asq))
        )
        text!(ax, minimum(kx_list) + label_shift, maximum(ky_list) - label_shift, text="(a)", color=:white, align=(:center, :center),)
        cross_mark(ax, k_1_max, 0.1)
        cross_mark(ax, k_2_max, 0.1)
        hidexdecorations!(ax)
        
        ax = Axis(f[1, 2], title=L"$\delta$ factor", titlefont=:regular, aspect=1)
        heatmap!(ax, kx_list, ky_list, E_residue, 
            colorrange=(minimum(E_residue), maximum(E_residue))
        )
        text!(ax, minimum(kx_list) + label_shift, maximum(ky_list) - label_shift, color=:white, align=(:center, :center),)
        cross_mark(ax, k_1_max, 0.1)
        cross_mark(ax, k_2_max, 0.1)
        hidexdecorations!(ax)
        hideydecorations!(ax)
    
        ax = Axis(f[1, 3], 
            # title=rich(
            #     "Final contribution from one ",
            #     rich("k", font=:italic),
            #     subscript("1")
            # ),
        aspect=1)
        Asq_with_constraint = E_residue .* Asq
        heatmap!(ax, kx_list, ky_list, Asq_with_constraint)
        text!(ax, minimum(kx_list) + label_shift, maximum(ky_list) - label_shift, color=:white, align=(:center, :center),)
        cross_mark(ax, k_1_max, 0.1)
        cross_mark(ax, k_2_max, 0.1)
        hidexdecorations!(ax)
        hideydecorations!(ax)
    
        #endregion
        ##########################################
    end
    
    # The k_1 = k_2 case
    let 
        k_1_max, k_2_max, ω = let k_1 = (m_e / M * (P - w) - k) / 2
            momenta = momentum_calc_ehh(trion, P, k, k_1)
            k_h1 = momenta.k_h1
            k_h2 = momenta.k_h2
            ω = E_trion_ehh(trion, P) - E_residue_ehh(trion, k_h1, k_h2)
            k_2 = momenta.k_2
            k_1, k_2, ω
        end
    
        ##########################################
        #region Calculating the k_1 = 0 case
        
        E_residue = map(k1_grid) do k_1
            k_h1 = momentum_calc_ehh(trion, P, k, k_1).k_h1
            k_h2 = momentum_calc_ehh(trion, P, k, k_1).k_h2 
            broaden(ω - E_trion_ehh(trion, P) + E_residue_ehh(trion, k_h1, k_h2))
        end
        
        Asq = map(k1_grid) do k_1
            k_2 = momentum_calc_ehh(trion, P, k, k_1).k_2
            abs2(A_k1k2(k_1, k_2))
        end
        
        #endregion
        ##########################################
        
        ##########################################
        #region Plotting the k_1 = 0 case
        
        ax = Axis(f[2, 1], #=title="Wave function structural factor",=# aspect=1)
        heatmap!(ax, kx_list, ky_list, Asq,
            colorrange=(minimum(Asq), maximum(Asq))
        )
        text!(ax, minimum(kx_list) + label_shift, maximum(ky_list) - label_shift, text="(b)", color=:white, align=(:center, :center),)
        cross_mark(ax, k_1_max, 0.1)
        cross_mark(ax, k_2_max, 0.1)
        hidexdecorations!(ax)
        #hideydecorations!(ax)
        
        ax = Axis(f[2, 2], #=title="Dispersion relation constraint",=# aspect=1)
        heatmap!(ax, kx_list, ky_list, E_residue, 
            colorrange=(minimum(E_residue), maximum(E_residue))
        )
        text!(ax, minimum(kx_list) + label_shift, maximum(ky_list) - label_shift, color=:white, align=(:center, :center),)
        cross_mark(ax, k_1_max, 0.1)
        cross_mark(ax, k_2_max, 0.1)
        hidexdecorations!(ax)
        hideydecorations!(ax)
    
        ax = Axis(f[2, 3], 
            # title=rich(
            #     "Final contribution from one ",
            #     rich("k", font=:italic),
            #     subscript("1")
            # ),
        aspect=1)
        Asq_with_constraint = E_residue .* Asq
        heatmap!(ax, kx_list, ky_list, Asq_with_constraint)
        text!(ax, minimum(kx_list) + label_shift, maximum(ky_list) - label_shift, color=:white, align=(:center, :center),)
        cross_mark(ax, k_1_max, 0.1)
        cross_mark(ax, k_2_max, 0.1)
        hidexdecorations!(ax)
        hideydecorations!(ax)
    
        #endregion
        ##########################################
    end
    
end

let 
    trion = Intervalley2DChandraTrion(
        m_h = m_h,
        m_e = m_e,
        w = w,
        E_g = E_g,
        E_B = 0.76,
        a = 10.3,
        b = 25.2
    )
    A_k1k2 = wfn(trion)
    
    P = 0.8w + SA[0.0, 0.1]

    let 
        k = SVector{2, Float64}([1.3, 0.0])
        k1_max, ω = let mom_report = momentum_calc_eeh_e2(trion, P, k, SA[0.0, 0.0])
            k_1 = mom_report.k_1
            k_e1 = mom_report.k_e1
            k_h = mom_report.k_h
            ω = E_trion_eeh(trion, P) - E_residue_eeh_e2(trion, k_e1, k_h)
            k_1, ω
        end
        
        E_residue_k1 = map(k1_grid) do k_1
            k_e1 = momentum_calc_eeh_e2(trion, P, k, k_1).k_e1
            k_h = momentum_calc_eeh_e2(trion, P, k, k_1).k_h
            broaden(ω - E_trion_eeh(trion, P) + E_residue_eeh_e2(trion, k_e1, k_h)) 
        end

        Asq_k1 = map(k1_grid) do (k_x, k_y)
            k_1 = SVector{2, Float64}([k_x, k_y])
            k_2 = momentum_calc_eeh_e2(trion, P, k, k_1).k_2
            abs2(A_k1k2(k_1, k_2))
        end

        aspect_ratio = (maximum(kx_list) - minimum(kx_list)) / (maximum(ky_list) - minimum(ky_list))
        ax = Axis(f[3, 1], aspect=aspect_ratio,)
        heatmap!(ax, kx_list, ky_list, Asq_k1)
        text!(ax, minimum(kx_list) + label_shift, maximum(ky_list) - label_shift, text="(c)", color=:white, align=(:center, :center),)
        hidexdecorations!(ax)
        #hideydecorations!(ax)
        
        cross_mark(ax, k1_max, 0.1)

        ax = Axis(f[3, 2], aspect=1)
        heatmap!(ax, kx_list, ky_list, E_residue_k1, colorrange=(minimum(E_residue_k1), maximum(E_residue_k1)))
        cross_mark(ax, k1_max, 0.1)
        hidexdecorations!(ax)
        hideydecorations!(ax)

        ax = Axis(f[3, 3], aspect=1)
        Asq_with_constraint = Asq_k1 .* E_residue_k1
        heatmap!(ax, kx_list, ky_list, Asq_with_constraint, colorrange=(minimum(Asq_with_constraint), maximum(Asq_with_constraint)))
        cross_mark(ax, k1_max, 0.1)
        hidexdecorations!(ax)
        hideydecorations!(ax)
    end

    let 
        M = (2m_e + m_h)
        k = (m_e + m_h) / M * w + m_e / M * P
        k1_max = let mom_report = momentum_calc_eeh_e2(trion, P, k, SA[0.0, 0.0])
            k_1 = mom_report.k_1
            k_e1 = mom_report.k_e1
            k_h = mom_report.k_h
            k_1
        end
        
        ω = 1.5
        
        E_residue_k1 = map(k1_grid) do k_1
            k_e1 = momentum_calc_eeh_e2(trion, P, k, k_1).k_e1
            k_h = momentum_calc_eeh_e2(trion, P, k, k_1).k_h
            broaden(ω - E_trion_eeh(trion, P) + E_residue_eeh_e2(trion, k_e1, k_h)) 
        end

        Asq_k1 = map(k1_grid) do (k_x, k_y)
            k_1 = SVector{2, Float64}([k_x, k_y])
            k_2 = momentum_calc_eeh_e2(trion, P, k, k_1).k_2
            abs2(A_k1k2(k_1, k_2))
        end

        aspect_ratio = (maximum(kx_list) - minimum(kx_list)) / (maximum(ky_list) - minimum(ky_list))
        ax = Axis(f[4, 1], aspect=aspect_ratio,)
        heatmap!(ax, kx_list, ky_list, Asq_k1)
        text!(ax, minimum(kx_list) + label_shift, maximum(ky_list) - label_shift, text="(d)", color=:white, align=(:center, :center),)
        #hidexdecorations!(ax)
        #hideydecorations!(ax)
        
        cross_mark(ax, k1_max, 0.1)

        ax = Axis(f[4, 2], aspect=1)
        heatmap!(ax, kx_list, ky_list, E_residue_k1, colorrange=(minimum(E_residue_k1), maximum(E_residue_k1)))
        cross_mark(ax, k1_max, 0.1)
        #hidexdecorations!(ax)
        hideydecorations!(ax)

        ax = Axis(f[4, 3], aspect=1)
        Asq_with_constraint = Asq_k1 .* E_residue_k1
        heatmap!(ax, kx_list, ky_list, Asq_with_constraint, colorrange=(minimum(Asq_with_constraint), maximum(Asq_with_constraint)))
        cross_mark(ax, k1_max, 0.1)
        #hidexdecorations!(ax)
        hideydecorations!(ax)
        
    end

    
end

colsize!(f.layout, 1, Relative(0.3))
colsize!(f.layout, 2, Relative(0.3))
colsize!(f.layout, 3, Relative(0.3))

Label(f[4, 1:3, Bottom()], rich(rich("k", subscript("x", font=:regular), font=:italic), " (Å⁻¹)"), tellwidth=false, tellheight=false, padding=(0, 0, 0, 5))
Label(f[1:4, 1, Left()], rich(rich("k", subscript("y", font=:regular), font=:italic), " (Å⁻¹)"), tellwidth=false, tellheight=false, padding=(0, 70, 0, 0), rotation=π/2)

colgap!(f.layout, 5)
rowgap!(f.layout, -60)

save("ehh-momentum-conservation-factor-compact-2.png", f)
f