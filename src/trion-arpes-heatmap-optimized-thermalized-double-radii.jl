using ProgressMeter
include("trion-solver.jl")
include("arpes.jl")

# Finite momentum trion, slow version 
@time let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_h = 0.8184,
    m_e = 0.4794,
    E_g = 0.7830,
    ϵ = 6.4, 
    E_B = -0.1,
    w = 0.7,
    β = 10,
    a = 10.3,
    b = 25.2,
    kx_points = LinRange(-1.2, 1.5, 100),
    ky_points = LinRange(-1.2, 1.5, 100),
    dkx = step(kx_points), 
    dky = step(ky_points),
    dk = dkx * dky,
    k_points = map(t -> SVector{2, Float64}(collect(t)), 
        collect(Iterators.product(kx_points, ky_points))
    ),
    ω_points = LinRange(0, 1.5, 500),
    P_Tx = 0.5,
    Qx_points =  LinRange(w - 0.5, w + 0.5, 20),
    Qy_points =  LinRange(- 0.5, 0.5, 20),
    Q_points = map(t -> SVector{2, Float64}(collect(t)), 
        collect(Iterators.product(Qx_points, Qy_points))
    ) 

    ϵ_v(k) =  - norm(k)^2 / 2m_h * inv_eV
    ϵ_c(k) =    norm(k - w)^2 / 2m_e * inv_eV + E_g
    
    M = 2m_h + m_e

    ham = IndirectTwoBandModel2D(m_e, m_h, E_g, SVector{2, Float64}([w, 0.0]))
    dielectric = Dielectric2D(ϵ)
    broadening(x) = @fastmath exp(- σ^2 * x^2)

    Z_partition = sum(Q_points) do Q_point
        exp(-β * (Q_point - [w, 0.0])' * (Q_point - [w, 0.0]) / 2M)
    end
    progress = Progress(length(Q_points)) 
    A_kω_β = sum(Q_points) do Q_point
        next!(progress)
        A_kω_Q = single_trion_arpes_signature_thread(
            ham, E_B, 
            ground_state_1s_double_radii(IndirectTwoBandTrion2D(ham, dielectric), a, b), 
            Q_point, 
            k_points, ω_points, 
            map(kx -> SVector{2, Float64}([kx, 0.0]), 
                kx_points
            ),
            broadening
        )   
        exp(-β * (Q_point - [w, 0.0])' * (Q_point - [w, 0.0]) / 2M) * A_kω_Q / Z_partition
    end  
    
    p = heatmap(kx_points, ω_points, A_kω_β', 
        c = arpes_colormap(transparency_gradience),
        colorbar = false, 
        legend_foreground_color = :transparent,
        legend = :topleft,
        dpi = 500, 
        grid = false, 
        framestyle = :box,
    )

    plot!(p, kx_points, ϵ_v.(kx_points), 
        label = "",
        linestyle = :dash,
        c = colorant"lightskyblue3")
    plot!(p, kx_points, ϵ_c.(kx_points), 
        label = "",
        linestyle = :dash,
        c = colorant"lightskyblue3")
        
    kx_near_valley = LinRange(-0.2, 1.2, 100)
    k_points_near_valley = map(x -> [x, 0.0], kx_near_valley)
    # Position of intensity peak: the brightest spot in one signature, 
    # with varying P_T 
    plot!(p, kx_near_valley, 
        inv_eV * norm.(k_points_near_valley .- [[w, 0.0]]).^2 * (M / (2 * m_e^2)) 
        .- inv_eV * norm.(k_points_near_valley .- [[w, 0.0]]).^2 * (m_h / m_e^2)
        .+ E_B .+ E_g,
        #label = raw"Varying $\mathbf{P}_\mathrm{T}$, $\mathbf{k}_1 = \mathbf{k}_2 = 0$",
        label = raw"Varying $\mathbf{P}$, $\mathbf{k}_{\mathrm{1}} = \mathbf{k}_{\mathrm{2}} = 0$",
        linestyle = :dot,
        c = colorant"firebrick2")
    # Dispersion relation when P_T is fixed to a constant 
    plot!(p, kx_near_valley, 
        ϵ_v.([[w, 0.0]] .- k_points_near_valley) 
        .+ E_B .+ E_g,
        #label = raw"fixed $\mathbf{P}_\mathrm{T}$, $\mathbf{k}_1 = 0$ or $\mathbf{k}_2 = 0$",
        label = raw"Fixed $\mathbf{P}$, $\mathbf{k}_{\mathrm{1}}$ or $\mathbf{k}_{\mathrm{2}}$ is $0$",
        linestyle = :dot, 
        c = colorant"aqua")
    plot!(p, kx_near_valley, 
        2 * ϵ_v.(([[w, 0.0]] .- k_points_near_valley) ./ 2) 
        .+ E_B .+ E_g,
        label = raw"Fixed $\mathbf{P}$, $\mathbf{k}_{\mathrm{1}} = \mathbf{k}_{\mathrm{2}} $",
        linestyle = :dot, 
        c = colorant"deepskyblue2")
    
    ylims!(p, (-0.8, 1.8))
    xlims!(p, (-0.8, 1.2))
    xlabel!(p, raw"$k$ (Å)")
    ylabel!(p, raw"$\omega$ (eV)")

    cb = heatmap(
        [1 0; 0 1], 
        clims=(0,1), 
        framestyle=:none, 
        c=arpes_colormap(white_gradience), 
        cbar=true, 
        lims=(-1,0)
    )

    p = custom_colorbar(p, cb)
    savefig(p, "trion-β-$β-w-$w.png")
end