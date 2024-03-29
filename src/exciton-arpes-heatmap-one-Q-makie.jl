using ProgressMeter
include("exciton-solver.jl")
include("arpes-makie.jl")

# Interestingly, when I average over the 2D grid, 
# the resulting signature loos less like the signature in Photoemission signature of excitons
# I don't know: is the latter obtained by confining the problem on a 1D line?
p = let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_v = 0.21,
    m_c = 0.37,
    E_g = 2.84,
    ϵ = 6.4, 
    E_B = -0.1,
    w = 0,#4π / (3 * 3.144817974),
    β = 10,
    broadening(x) = exp(- σ^2 * x^2),
    kx_points = LinRange(-1.2, 1.5, 100),
    ky_points = LinRange(-1.2, 1.5, 100),
    k_points = map(t -> SVector{2, Float64}(collect(t)), 
        collect(Iterators.product(kx_points, ky_points))
    ),
    ω_points = LinRange(0, 3, 500),
    hexagonal_edge = 4π / (3 * 3.144817974),
    Q = 0.2 * hexagonal_edge,
    Q_point = SVector{2, Float64}([Q, 0.0])

    ϵ_v(k) = - norm(k)^2 / 2m_v * inv_eV
    ϵ_c(k) = E_g + norm(k - w)^2 / 2m_c * inv_eV

    M = m_c + m_v
    α = m_c / M
    μ = m_c * m_v / (m_c + m_v)
    E_SQ(Q) = E_g + E_B + (Q - w)^2/2M * inv_eV

    ham = IndirectTwoBandModel2D(m_c, m_v, E_g, SVector{2, Float64}([w, 0.0]))
    dielectric = Dielectric2D(ϵ)

    A_kω_Q = single_exciton_arpes_signature_def(
        ham, E_B, 
        ground_state_1s(IndirectTwoBandExciton2D(ham, dielectric)), 
        Q_point, 
        ω_points, 
        map(kx -> SVector{2, Float64}([kx, 0.0]), 
            kx_points
        ),
        broadening
    )   
    
    fig = Figure()
    ax_heatmap = Axis(fig[1, 1], 
        xlabel = L"$k$ (Å)", 
        ylabel = L"$ω$ (eV)",
        xlabelsize = 20,
        ylabelsize = 20,
        xticklabelsize = 18,
        yticklabelsize = 18,
    )

    heatmap!(ax_heatmap, kx_points, ω_points, A_kω_Q, 
        colormap = arpes_colormap(transparency_gradience),
    ) 
    lines!(ax_heatmap, kx_points, ϵ_v.(kx_points), 
        linestyle = :dash,
        colormap = [colorant"lightskyblue3"],
        color = 1,
        colorrange = (0, 10),
    )
    lines!(ax_heatmap, kx_points, ϵ_c.(kx_points), 
        linestyle = :dash,
        colormap = [colorant"lightskyblue3"],
        color = 1,
        colorrange = (0, 10),
    )

    k_points_near_valley = LinRange(-0.2, 1.2, 100)
    lines!(ax_heatmap, k_points_near_valley, ϵ_c.(k_points_near_valley) .+ E_B,
        label = L"Varying $\mathbf{Q}$",
        linestyle = :dot,
        colormap = [colorant"firebrick2"],
        color = 1,
        colorrange = (0, 10),
    )
    lines!(ax_heatmap, k_points_near_valley, ϵ_v.(k_points_near_valley .- Q_point[1]) .+ E_SQ(Q),
        label = L"Fixed $\mathbf{Q}$",
        linestyle = :dot,
        colormap = [colorant"aqua"],
        color = 1,
        colorrange = (0, 10),
    )
    
    ylims!(ax_heatmap, (-0.8, 4))
    xlims!(ax_heatmap, (-0.8, 1.2)) 
    hidedecorations!(ax_heatmap, ticks = false, ticklabels = false, label = false)
    
    axislegend(ax_heatmap, backgroundcolor = RGBA(1, 1, 1, 0.5), framecolor = RGBA(1, 1, 1, 0.5), labelsize = 18)
    Colorbar(fig[1, 2], 
        colormap = arpes_colormap(transparency_gradience),
        colorrange = (minimum(A_kω_Q), maximum(A_kω_Q))
    )
    
    save("exciton-Q-$(round(Q / hexagonal_edge, digits = 3))-w-$(round(w / hexagonal_edge, digits = 3))-makie.png", fig)
    save("exciton-Q-$(round(Q / hexagonal_edge, digits = 3))-w-$(round(w / hexagonal_edge, digits = 3))-makie.pdf", fig)
    fig
end
