# Plot the contribution of the trion wave function structural factor and the energy conservation delta function

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

function cross_mark(ax::Axis, pos::AbstractVector, len::Float64)
    lines!(ax, [pos[1] - len, pos[1] + len], [pos[2], pos[2]], color=:white, linewidth=0.5)
    lines!(ax, [pos[1], pos[1]], [pos[2] - len, pos[2] + len], color=:white, linewidth=0.5)
end

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
M = 2m_h + m_e
A_k1k2 = wfn(trion)

k = SVector{2, Float64}([0.3, -0.1])
P = 1.2w + SA[0.0, 0.1]

kx_list = LinRange(-0.5, 0.5, 200)
ky_list = LinRange(-0.5, 0.5, 200)
k1_grid = map(Iterators.product(kx_list, ky_list)) do (k_x, k_y)
    SA[k_x, k_y]
end
ω_list = LinRange(0, 3.0, 200) #LinRange(-8, 5, 200)

set_theme!(fontsize=20)
f = Figure(size=(700, 1000))

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
    
    ax = Axis(f[1, 1], #=title="Wave function structural factor",=# aspect=1)
    Label(f[1, 1, TopLeft()], "(a)", padding=(0, 15, 15, 0))
    heatmap!(ax, kx_list, ky_list, Asq,
        colorrange=(minimum(Asq), maximum(Asq))
    )
    cross_mark(ax, k_1_max, 0.1)
    cross_mark(ax, k_2_max, 0.1)
    
    ax = Axis(f[2, 1], #=title="Dispersion relation constraint",=# aspect=1)
    Label(f[2, 1, TopLeft()], "(b)", padding=(0, 15, 15, 0))
    heatmap!(ax, kx_list, ky_list, E_residue, 
        colorrange=(minimum(E_residue), maximum(E_residue))
    )
    cross_mark(ax, k_1_max, 0.1)
    cross_mark(ax, k_2_max, 0.1)

    ax = Axis(f[3, 1], 
        # title=rich(
        #     "Final contribution from one ",
        #     rich("k", font=:italic),
        #     subscript("1")
        # ),
    aspect=1)
    Label(f[3, 1, TopLeft()], "(c)", padding=(0, 15, 15, 0))
    Asq_with_constraint = E_residue .* Asq
    heatmap!(ax, kx_list, ky_list, Asq_with_constraint)
    cross_mark(ax, k_1_max, 0.1)
    cross_mark(ax, k_2_max, 0.1)

    #endregion
    ##########################################

    # The width of column 1 should be the same as the height of row 1
    colsize!(f.layout, 1, Aspect(1, 1))
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
    
    ax = Axis(f[1, 2], #=title="Wave function structural factor",=# aspect=1)
    Label(f[1, 2, TopLeft()], "(d)", padding=(0, 15, 15, 0))
    heatmap!(ax, kx_list, ky_list, Asq,
        colorrange=(minimum(Asq), maximum(Asq))
    )
    cross_mark(ax, k_1_max, 0.1)
    cross_mark(ax, k_2_max, 0.1)
    
    ax = Axis(f[2, 2], #=title="Dispersion relation constraint",=# aspect=1)
    Label(f[2, 2, TopLeft()], "(e)", padding=(0, 15, 15, 0))
    heatmap!(ax, kx_list, ky_list, E_residue, 
        colorrange=(minimum(E_residue), maximum(E_residue))
    )
    cross_mark(ax, k_1_max, 0.1)
    cross_mark(ax, k_2_max, 0.1)

    ax = Axis(f[3, 2], 
        # title=rich(
        #     "Final contribution from one ",
        #     rich("k", font=:italic),
        #     subscript("1")
        # ),
    aspect=1)
    Label(f[3, 2, TopLeft()], "(f)", padding=(0, 15, 15, 0))
    Asq_with_constraint = E_residue .* Asq
    heatmap!(ax, kx_list, ky_list, Asq_with_constraint)
    cross_mark(ax, k_1_max, 0.1)
    cross_mark(ax, k_2_max, 0.1)

    #endregion
    ##########################################

    # The width of column 2 should be the same as the height of row 1
    colsize!(f.layout, 2, Aspect(1, 1))
end

save("momentum-conservation-factor.png", f)
f