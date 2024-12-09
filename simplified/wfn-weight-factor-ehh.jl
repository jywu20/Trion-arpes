using CairoMakie
include("wfn.jl")
include("ehh.jl")
include("broadening.jl")
include("units.jl")

let f = Figure()
    
    kx_list = LinRange(-0.5, 1.3, 100)
    ky_list = LinRange(-0.5, 0.5, 100)
    #kx_list = LinRange(-1, 1, 100)
    #ky_list = LinRange(-1, 1, 100)

    ax = Axis(f[1, 1], aspect=maximum(kx_list) / maximum(ky_list))

    k1_list = Iterators.product(kx_list, ky_list)
    w_side = 4π / (3 * 3.144817974)
    w = SVector{2, Float64}([w_side, 0.0])

    trion = Intervalley2DChandraTrion(
        m_h = 0.21,
        m_e = 0.37,
        w = w,
        E_g = 2.84,
        E_B = 0.1,
        a = 10.3,
        b = 25.2
    )
    A_k1k2 = wfn(trion)

    k = SVector{2, Float64}([-1, 0.0])
    P = w

    Asq = map(k1_list) do (k_x, k_y)
        k_1 = SVector{2, Float64}([k_x, k_y])
        k_2 = momentum_calc_ehh(trion, P, k, k_1).k_2
        abs2(A_k1k2(k_1, k_2))
    end
    
    heatmap!(ax, kx_list, ky_list, Asq)
    
    ω = 4.0
    σ = 20.0
    m_e = trion.m_e
    m_h = trion.m_h
    E_g = trion.E_g
    E_B = trion.E_B
    
    broaden = gaussian_broadening(σ)

    ax = Axis(f[2, 1], aspect=maximum(kx_list) / maximum(ky_list))
    E_residue = map(k1_list) do (k_x, k_y)
        k_1 = SVector{2, Float64}([k_x, k_y])
        k_h1 = momentum_calc_ehh(trion, P, k, k_1).k_h1
        k_h2 = momentum_calc_ehh(trion, P, k, k_1).k_h2 
        broaden(ω - E_trion(trion, P) + E_residue_ehh(trion, k_h1, k_h2))
    end
    
    heatmap!(ax, kx_list, ky_list, E_residue)
    

    f
end

let f = Figure()
    ax = Axis(f[1, 1])
    ϕ(k_1, k_2) = ϕ_1s(SA[k_1, 0.0], 10.3) * ϕ_1s(SA[k_2, 0.0], 25.2) + ϕ_1s(SA[k_2, 0.0], 10.3) * ϕ_1s(SA[k_1, 0.0], 25.2)
    
    k_1_list = LinRange(-10, 10, 100)
    ϕ_list = map(k_1_list) do k_1
        k_2 = 5 - k_1
        ϕ(k_1, k_2)
    end
    lines!(ax, k_1_list, ϕ_list)

    f
end