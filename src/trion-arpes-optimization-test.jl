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
    kx_points = LinRange(-1.2, 1.2, 100),
    ky_points = LinRange(-1.2, 1.2, 100),
    dkx = step(kx_points), 
    dky = step(ky_points),
    dk = dkx * dky,
    k_points = map(t -> SVector{2, Float64}(collect(t)), 
        collect(Iterators.product(kx_points, ky_points))
    ),
    ω_points = LinRange(0, 1.5, 500),
    P_Tx = 0.8,
    Q_point  = SVector{2, Float64}([P_Tx, 0.0]) 

    ϵ_v(k) =  - norm(k)^2 / 2m_h * inv_eV
    ϵ_c(k) =    norm(k - w)^2 / 2m_e * inv_eV + E_g
    
    M = 2m_h + m_e

    ham = IndirectTwoBandModel2D(m_e, m_h, E_g, SVector{2, Float64}([w, 0.0]))
    dielectric = Dielectric2D(ϵ)
    broadening(x) = @fastmath exp(- σ^2 * x^2)

    @time A_kω_Q = single_trion_arpes_signature_thread(
        ham, E_B, 
        ground_state_1s(ham, dielectric), 
        Q_point, 
        k_points, ω_points, 
        map(kx -> SVector{2, Float64}([kx, 0.0]), 
            kx_points
        ),
        broadening
    )     
    
    @time A_kω_Q_def = single_trion_arpes_signature_def(
        ham, E_B, 
        ground_state_1s(ham, dielectric), 
        Q_point, 
        k_points, ω_points, 
        map(kx -> SVector{2, Float64}([kx, 0.0]), 
            kx_points
        ),
        broadening
    )     
    
    sum(abs.(A_kω_Q_def - A_kω_Q)) 
end