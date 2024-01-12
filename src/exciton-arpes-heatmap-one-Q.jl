using ProgressMeter
include("exciton-solver.jl")
include("arpes.jl")

# Interestingly, when I average over the 2D grid, 
# the resulting signature loos less like the signature in Photoemission signature of excitons
# I don't know: is the latter obtained by confining the problem on a 1D line?
p = let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_v = 0.8184,
    m_c = 0.4794,
    E_g = 0.7830,
    ϵ = 6.4, 
    E_B = -0.1,
    w = 0.7, # The value used in the literature: 0.7,
    β = 10,
    broadening(x) = exp(- σ^2 * x^2),
    kx_points = LinRange(-1.2, 1.5, 100),
    ky_points = LinRange(-1.2, 1.5, 100),
    k_points = map(t -> SVector{2, Float64}(collect(t)), 
        collect(Iterators.product(kx_points, ky_points))
    ),
    ω_points = LinRange(-1.2, 1.2, 100),
    Q = 0.7,
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
    
    p = heatmap(kx_points, ω_points, A_kω_Q', 
        c = arpes_colormap(transparency_gradience),
        colorbar = false,
        legend_foreground_color = :transparent,
        legend = :topleft,
        dpi = 500, 
        grid = false, 
        framestyle = :box,) 
    plot!(p, kx_points, ϵ_v.(kx_points), 
        label = "",
        linestyle = :dash,
        c = colorant"lightskyblue3")
    plot!(p, kx_points, ϵ_c.(kx_points), 
        label = "",
        linestyle = :dash,
        c = colorant"lightskyblue3")

    k_points_near_valley = LinRange(-0.2, 1.2, 100)
    plot!(p, k_points_near_valley, ϵ_c.(k_points_near_valley) .+ E_B,
        label = raw"Varying $\mathbf{Q}$",
        linestyle = :dot,
        c = colorant"firebrick2")
    plot!(p, k_points_near_valley, ϵ_v.(k_points_near_valley .- Q_point[1]) .+ E_SQ(Q),
        label = raw"Fixed $\mathbf{Q}$",
        linestyle = :dot,
        c = colorant"aqua")
    
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
    savefig(p, "exciton-Q-$Q-w-$w.png")
end
