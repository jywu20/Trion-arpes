using ProgressMeter
include("trion-solver.jl")
include("arpes-makie.jl")

@time let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_h = 0.21,
    m_e = 0.37,
    E_g = 2.84,
    ϵ = 6.4, 
    E_B = -0.1,
    w = 0.7,
    β = 1,
    a = 10.3,
    b = 25.2,
    kx_points = LinRange(-1.2, 1.2, 200),
    ky_points = LinRange(-1.2, 1.2, 200),
    dkx = step(kx_points), 
    dky = step(ky_points),
    dk = dkx * dky,
    k_points = map(t -> SVector{2, Float64}(collect(t)), 
        collect(Iterators.product(kx_points, ky_points))
    ),
    ω_points = LinRange(-1, 3, 500),
    P_x = 0.7,
    P_point  = SVector{2, Float64}([P_x, 0.0]) 

    ϵ_v2(k) =  - norm(k - w)^2 / 2m_h * inv_eV
    ϵ_v1(k) =  - norm(k)^2     / 2m_h * inv_eV
    ϵ_c1(k) =    norm(k)^2     / 2m_e * inv_eV + E_g
    ϵ_c2(k) =    norm(k - w)^2 / 2m_e * inv_eV + E_g
    E_SP(P) = inv_eV * norm(P .- w)^2 / 2M .+ E_B .+ 2E_g
    
    M = 2m_h + m_e

    ham = IndirectTwoBandModel2D(m_e, m_h, E_g, SVector{2, Float64}([w, 0.0]))
    dielectric = Dielectric2D(ϵ)
    broadening(x) = @fastmath exp(- σ^2 * x^2)

    A_kω_Q = trionarpes_e1e2h1_thread(
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
    lines!(ax_heatmap, kx_points, ϵ_c1.(kx_points), 
        linestyle = :dash,
        colormap = [colorant"lightskyblue3"], 
        color = 1, 
        colorrange = (0, 10),
    )
    lines!(ax_heatmap, kx_points, ϵ_c2.(kx_points), 
        linestyle = :dash,
        colormap = [colorant"lightskyblue3"], 
        color = 1, 
        colorrange = (0, 10),
    )
        
    kx_near_valley = LinRange(-0.7, 0.7, 100)
    k_points_near_valley = map(x -> [x, 0.0], kx_near_valley)
    # k_1 when e1 is out  
    k1_near_valley_e1 = kx_near_valley .+ m_e / M * (w - P_x) 
    k2_near_valley_e1 = - m_e / (m_e + m_h) * k1_near_valley_e1
    k2_near_valley_e2 = kx_near_valley .- m_e / M * P_x .- (m_e + m_h) / M * w
    k1_near_valley_e2 = - m_e / (m_e + m_h) * k2_near_valley_e2
    # K_1 when k_2 = 0
    k1_near_valley_extreme = (P_x - w) * m_e / M .- kx_near_valley
    # When k1, k2 = 0, the values of kh1 and kh2 w.r.t. Q 
    Q_wrt_ke = w .+ M / m_e * kx_near_valley
    kh1_k1k2max = m_h / M * Q_wrt_ke .- m_h / M * w 
    kh2_k1k2max = m_h / M * Q_wrt_ke .+ (m_e + m_h) / M * w

    # Position of intensity peak: the brightest spot in one signature, 
    # with varying P_T 
#    lines!(p, kx_near_valley, 
#        ϵ_v1.(kh1_k1k2max) + ϵ_v2.(kh2_k1k2max) .+ E_SP.(Q_wrt_ke),
#        #label = raw"Varying $\mathbf{P}_\mathrm{T}$, $\mathbf{k}_1 = \mathbf{k}_2 = 0$",
#        label = raw"Varying $\mathbf{P}$, $\mathbf{k}_{\mathrm{1}} = \mathbf{k}_{\mathrm{2}} = 0$",
#        linestyle = :dot,
#        c = colorant"firebrick2")
        
    # Dispersion relation when P_T is fixed to a constant 
    #
    # e1 is driven out, leaving e2 and h; 
    # the total kinetic energy is minimized 
    lines!(ax_heatmap, kx_near_valley, 
        - ϵ_c2.(k2_near_valley_e1 .+ (m_e + m_h) / M * w .+ m_e / M * P_x) 
        .+ ϵ_v1.(m_h / M * (P_x - w) .- k1_near_valley_e1 - k2_near_valley_e1)
        .+ E_SP(P_x), 
        label = L"Fixed $\mathbf{P}$, $\mathrm{min} E_{\mathrm{kin}}$",
        linestyle = :dot, 
        colormap = [colorant"aqua"], 
        color = 1, 
        colorrange = (0, 10),
    )
    # e1 is driven out, leaving e2 and h; 
    # k2 is set to zero.
    lines!(ax_heatmap, kx_near_valley, 
        - ϵ_c2.((m_e + m_h) / M * w .+ m_e / M * P_x) 
        .+ ϵ_v1.(m_h / M * (P_x - w) .- k1_near_valley_e1)
        .+ E_SP(P_x), 
        label = L"Fixed $\mathbf{P}$, $\mathbf{k}_{\mathrm{1}}$ or $\mathbf{k}_{\mathrm{2}}$ is $0$",
        linestyle = :dot, 
        colormap = [colorant"deepskyblue2"], 
        color = 1, 
        colorrange = (0, 10),
    )
    # e2 is driven out, leaving e1 and h; 
    # the total kinetic energy is minimized 
    lines!(ax_heatmap, kx_near_valley, 
        - ϵ_c1.(k1_near_valley_e2 .- m_e / M * w .+ m_e / M * P_x) 
        .+ ϵ_v1.(m_h / M * (P_x - w) .- k2_near_valley_e2 .- k1_near_valley_e2)
        .+ E_SP(P_x), 
        label = L"Fixed $\mathbf{P}$, $\mathrm{min} E_{\mathrm{kin}}$",
        linestyle = :dot, 
        colormap = [colorant"aqua"], 
        color = 1, 
        colorrange = (0, 10),
    )
    # e2 is driven out, leaving e1 and h; 
    # k1 is set to zero 
    lines!(ax_heatmap, kx_near_valley, 
        - ϵ_c1.( -m_e / M * w .+ m_e / M * P_x) 
        .+ ϵ_v1.(m_h / M * (P_x - w) .- k2_near_valley_e2)
        .+ E_SP(P_x), 
        label = L"Fixed $\mathbf{P}$, $\mathbf{k}_{\mathrm{1}}$ or $\mathbf{k}_{\mathrm{2}}$ is $0$",
        linestyle = :dot, 
        colormap = [colorant"deepskyblue2"], 
        color = 1, 
        colorrange = (0, 10),
    )

    axislegend(ax_heatmap, backgroundcolor = RGBA(1, 1, 1, 0.5), framecolor = RGBA(1, 1, 1, 0.5), labelsize = 18)

    ylims!(ax_heatmap, (-0.8, 4))
    xlims!(ax_heatmap, (-0.8, 1.2))
    hidedecorations!(ax_heatmap, ticks = false, ticklabels = false, label = false)
    
    save("trion-e1e2h1-P-$(P_point[1])-w-$w-makie.png", fig)
    fig
end

# For debugging
@time let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_h = 0.21,
    m_e = 0.37,
    E_g = 2.84,
    ϵ = 6.4, 
    E_B = -0.1,
    w = 0.7,
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
    ω_points = LinRange(-1, 1.5, 500),
    P_x = 0.7,
    P_point  = SVector{2, Float64}([P_x, 0.0]) 

    ϵ_v2(k) =  - norm(k - w)^2 / 2m_h * inv_eV
    ϵ_v1(k) =  - norm(k)^2     / 2m_h * inv_eV
    ϵ_c1(k) =    norm(k)^2     / 2m_e * inv_eV + E_g
    ϵ_c2(k) =    norm(k - w)^2 / 2m_e * inv_eV + E_g
    E_SP(P) = inv_eV * norm(P .- w)^2 / 2M .+ E_B .+ 2E_g
    
    M = 2m_h + m_e

    ham = IndirectTwoBandModel2D(m_e, m_h, E_g, SVector{2, Float64}([w, 0.0]))
    dielectric = Dielectric2D(ϵ)
    broadening(x) = @fastmath exp(- σ^2 * x^2)

    ϕ_ex = phi1sa1sb(IndirectTwoBandMat2D(ham, dielectric), a, b)
        
    kx_near_valley = LinRange(-0.7, 0.7, 100)
    k_points_near_valley = map(x -> [x, 0.0], kx_near_valley)
    # k_1 when e1 is out  
    k1_near_valley_e1 = kx_near_valley .+ m_e / M * (w - P_x) 
    k2_near_valley_e2 = kx_near_valley .- m_e / M * P_x .- (m_e + m_h) / M * w
    # K_1 when k_2 = 0
    k1_near_valley_extreme = (P_x - w) * m_e / M .- kx_near_valley
    # When k1, k2 = 0, the values of kh1 and kh2 w.r.t. Q 
    Q_wrt_ke = w .+ M / m_e * kx_near_valley
    kh1_k1k2max = m_h / M * Q_wrt_ke .- m_h / M * w 
    kh2_k1k2max = m_h / M * Q_wrt_ke .+ (m_e + m_h) / M * w

    # e1 is driven out, leaving e2 and h; 
    # k2 is set to zero.
    ke2 = k1_near_valley_e1 .+ (m_e + m_h) / M * w .+ m_e / M * P_x
    kh  = m_h / M * (P_x - w) .- 2k1_near_valley_e1
    plot(kx_near_valley, ϕ_ex.(ke2, kh))
end