using ProgressMeter
include("trion-solver.jl")
include("arpes-makie.jl")

# Finite momentum trion, slow version 
@time let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_h = 0.21,
    m_e = 0.37,
    E_g = 2.84,
    ϵ = 6.4, 
    E_B = -0.1,
    hexagonal_edge = 4π / (3 * 3.144817974),
    w = hexagonal_edge,
    β = 1,
    a = 10.3,
    b = 25.2,
    kx_points = LinRange(-1.2, 1.2, 100),
    ky_points = LinRange(-1.2, 1.2, 100),
    dkx = step(kx_points), 
    dky = step(ky_points),
    dk = dkx * dky,
    k_points = map(t -> SVector{2, Float64}(collect(t)), 
        collect(Iterators.product(kx_points, ky_points))
    ),
    ω_points = LinRange(0, 3, 500),
    P_Tx = 1.2 * w,
    P_point  = SVector{2, Float64}([P_Tx, 0.0]) 

    ϵ_v2(k) =  - norm(k - w)^2 / 2m_h * inv_eV
    ϵ_v1(k) =  - norm(k)^2     / 2m_h * inv_eV
    ϵ_c(k)  =    norm(k)^2     / 2m_e * inv_eV + E_g
    E_SP(P) = inv_eV * norm(P .- w)^2 / 2M .+ E_B .+ E_g
    
    M = 2m_h + m_e

    ham = IndirectTwoBandModel2D(m_e, m_h, E_g, SVector{2, Float64}([w, 0.0]))
    dielectric = Dielectric2D(ϵ)
    broadening(x) = @fastmath exp(- σ^2 * x^2)

    A_kω_Q = trionarpes_e1h1h2_thread(
        ham, E_B, 
        phi1sa1sb(IndirectTwoBandMat2D(ham, dielectric), a, b), 
        P_point, 
        k_points, ω_points, 
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
    
    lines!(ax_heatmap, kx_points, ϵ_v1.(kx_points), 
        linestyle = :dash,
        colormap = [colorant"lightskyblue3"], 
        color = 1, 
        colorrange = (0, 10),
    )
    lines!(ax_heatmap, kx_points, ϵ_v2.(kx_points), 
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
        
    kx_near_valley = LinRange(-0.7, 0.7, 100)
    k_points_near_valley = map(x -> [x, 0.0], kx_near_valley)
    # k_1 when k_1 = k_2
    k1_near_valley_equal   = 1/2 * ((P_Tx - w) * m_e / M .- kx_near_valley)
    # K_1 when k_2 = 0
    k1_near_valley_extreme = (P_Tx - w) * m_e / M .- kx_near_valley
    # When k1, k2 = 0, the values of kh1 and kh2 w.r.t. Q 
    Q_wrt_ke = w .+ M / m_e * kx_near_valley
    kh1_k1k2max = m_h / M * Q_wrt_ke .- m_h / M * w 
    kh2_k1k2max = m_h / M * Q_wrt_ke .+ (m_e + m_h) / M * w

    # Position of intensity peak: the brightest spot in one signature, 
    # with varying P_T 
    lines!(ax_heatmap, kx_near_valley, 
        ϵ_v1.(kh1_k1k2max) + ϵ_v2.(kh2_k1k2max) .+ E_SP.(Q_wrt_ke),
        #label = raw"Varying $\mathbf{P}_\mathrm{T}$, $\mathbf{k}_1 = \mathbf{k}_2 = 0$",
        label = L"Varying $\mathbf{P}$, $\mathbf{k}_{\mathrm{1}} = \mathbf{k}_{\mathrm{2}} = 0$",
        linestyle = :dot,
        colormap = [colorant"firebrick2"],
        color = 1,
        colorrange = (0, 10),
    )
    println(ϵ_v1.(kh1_k1k2max) + ϵ_v2.(kh2_k1k2max) .+ E_SP(Q_wrt_ke))
    # Dispersion relation when P_T is fixed to a constant 
    lines!(ax_heatmap, kx_near_valley, 
        ϵ_v1.(k1_near_valley_extreme .- (w - P_Tx) * m_h / M) 
        .+ ϵ_v2((m_e + m_h) / M * w + m_h / M * P_Tx)
        .+ E_SP(P_Tx),
        #label = raw"fixed $\mathbf{P}_\mathrm{T}$, $\mathbf{k}_1 = 0$ or $\mathbf{k}_2 = 0$",
        label = L"Fixed $\mathbf{P}$, $\mathbf{k}_{\mathrm{1}}$ or $\mathbf{k}_{\mathrm{2}}$ is $0$",
        linestyle = :dot, 
        colormap = [colorant"aqua"],
        color = 1,
        colorrange = (0, 10),
    )
    lines!(ax_heatmap, kx_near_valley, 
        2 * ϵ_v1.(k1_near_valley_equal) 
        .+ E_SP(P_Tx),
        label = L"Fixed $\mathbf{P}$, $\mathbf{k}_{\mathrm{1}} = \mathbf{k}_{\mathrm{2}} $",
        linestyle = :dot, 
        colormap = [colorant"deepskyblue2"],
        color = 1,
        colorrange = (0, 10),
    )
    
    axislegend(ax_heatmap, backgroundcolor = RGBA(1, 1, 1, 0.5), framecolor = RGBA(1, 1, 1, 0.5), labelsize = 18)
    
    ylims!(ax_heatmap, (-0.8, 4))
    xlims!(ax_heatmap, (-0.8, 1.2))
    hidedecorations!(ax_heatmap, ticks = false, ticklabels = false, label = false)

    Colorbar(fig[1, 2], 
        colormap = arpes_colormap(transparency_gradience),
        colorrange = (minimum(A_kω_Q), maximum(A_kω_Q))
    )

    save("trion-e1h1h2-P-$(round(P_point[1] / hexagonal_edge, digits = 3))-w-$(round(w / hexagonal_edge, digits = 3))-makie.png", fig)
    save("trion-e1h1h2-P-$(round(P_point[1] / hexagonal_edge, digits = 3))-w-$(round(w / hexagonal_edge, digits = 3))-makie.pdf", fig)
    fig
end