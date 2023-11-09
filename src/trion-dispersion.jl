using Plots
using LinearAlgebra
using ProgressMeter
using Colors
using ColorSchemes
using LaTeXStrings
using StaticArrays

const fs = 1.5193
const eV = 3.80998

function transparency_gradience(color; steps = 100)
    α = alpha(color)
    r = red(color)
    b = blue(color)
    g = green(color)
    
    map(LinRange(0.0, α, steps)) do α_new
        RGBA(r, g, b, α_new)
    end
end

function white_gradience(color; steps = 100)
    ColorScheme([colorant"white", color])[LinRange(0, 1, steps)]
end

function arpes_colormap(gradient)
    rainbow = cgrad(:rainbow)
    starting_region = gradient(first(rainbow), steps = 50) 
    middle_region = rainbow[LinRange(0, 1, 500)]
    end_region = reverse(gradient(last(rainbow), steps = 400))
    
    PlotUtils.ContinuousColorGradient([starting_region..., middle_region..., end_region...])
end

function custom_colorbar(p, cb)
    l = @layout [a{0.95w} b]
    plot(p, cb, layout = l)
end

# The trion band index S and momentum Q should be pre-determined
function single_trion_arpes_signature(P_T, A_SQ_k1k2, m_c, m_v, E_g, E_B, broadening, k_points, ω_points, k_points_ARPES)
    m_T = 2m_v + m_c
    E_S_PT = eV * norm(P_T)^2 / 2m_T
    heatmap_points = Iterators.product(k_points_ARPES, ω_points) 

    A_kω = map(heatmap_points) do (k_e, ω)
        sum(k_points) do k_h1
            k_h2 = P_T - k_e - k_h1
            ϵ_v_kh1 = - eV * norm(k_h1)^2 / 2m_v
            ϵ_v_kh2 = - eV * norm(k_h2)^2 / 2m_v

            δ_factor = broadening(ω - E_S_PT - ϵ_v_kh1 - ϵ_v_kh2 - E_B - E_g)
            
            k_1 = k_h1 - m_v / m_T * P_T
            k_2 = k_h2 - m_v / m_T * P_T
            trion_structure_factor = abs(A_SQ_k1k2(k_1, k_2))^2
            
            δ_factor * trion_structure_factor
        end
    end
    
    A_kω
end


# Q = 0 trion 
let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_v = 0.8184,
    m_c = 0.4794,
    E_g = 0.7830,
    ϵ = 6.4, 
    E_B = -0.1,
    w = 0.7,
    β = 10,
    broadening(x) = exp(- σ^2 * x^2),
    kx_points = LinRange(-1.2, 1.2, 100),
    ky_points = LinRange(-1.2, 1.2, 100),
    dkx = step(kx_points), 
    dky = step(ky_points),
    dk = dkx * dky,
    k_points = map(collect, 
        collect(Iterators.product(kx_points, ky_points))
    ),
    ω_points = LinRange(0, 1.2, 500),
    Q_point  = SVector{2, Float64}([0.0, 0.0]) 

    ϵ_v(k) =  - norm(k)^2 / 2m_v * eV
    
    μ = m_c * m_v / (m_c + m_v)
    a_ex = 0.529 * ϵ / μ
    ϕ_1s(k) = 1 / (1 + norm(k)^2 * a_ex^2 / 4)^2
    A_SQ_k1k2(k_1, k_2) = ϕ_1s(k_1) * ϕ_1s(k_2)

    A_kω_Q = single_trion_arpes_signature(
        Q_point, A_SQ_k1k2, 
        m_c, m_v, E_g, E_B, broadening, k_points, ω_points, 
        map(kx -> SVector{2, Float64}([kx, 0.0]), kx_points))     
    
    p = heatmap(kx_points, ω_points, A_kω_Q', 
        c = arpes_colormap(transparency_gradience),
        colorbar = false) 
        
    k_points_near_valley = LinRange(-0.5, 0.5, 100)
    plot!(p, k_points_near_valley, ϵ_v.(k_points_near_valley) .+ E_B .+ E_g,
        legend = false,
        linestyle = :dot,
        c = colorant"aqua")
    # When Q = 0, and k_1 = k_h1 = - k_e / 2, k_2 = k_h2 = - k_e / 2
    plot!(p, k_points_near_valley, ϵ_v.(k_points_near_valley) / 2 .+ E_B .+ E_g,
        legend = false,
        linestyle = :dot,
        c = colorant"aqua")

    cb = heatmap(
        [1 0; 0 1], 
        clims=(0,1), 
        framestyle=:none, 
        c=arpes_colormap(white_gradience), 
        cbar=true, 
        lims=(-1,0)
    )
    custom_colorbar(p, cb)
end

# Finite momentum trion 
@time let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_v = 0.8184,
    m_c = 0.8, #0.4794,
    E_g = 0.7830,
    ϵ = 6.4, 
    E_B = -0.1,
    w = 0.7,
    β = 10,
    broadening(x) = exp(- σ^2 * x^2),
    kx_points = LinRange(-1.2, 1.2, 100),
    ky_points = LinRange(-1.2, 1.2, 100),
    dkx = step(kx_points), 
    dky = step(ky_points),
    dk = dkx * dky,
    k_points = map(t -> SVector{2, Float64}(collect(t)), 
        collect(Iterators.product(kx_points, ky_points))
    ),
    ω_points = LinRange(0, 1.5, 500),
    P_Tx = 0.0,
    Q_point  = SVector{2, Float64}([P_Tx, 0.0]) 

    ϵ_v(k) =  - norm(k)^2 / 2m_v * eV
    ϵ_c(k) =    norm(k)^2 / 2m_c * eV + E_g
    
    m_T = 2m_v + m_c

    μ = m_c * m_v / (m_c + m_v)
    a_ex = 0.529 * ϵ / μ
    ϕ_1s(k) = 1 / (1 + norm(k)^2 * a_ex^2 / 4)^2
    A_SQ_k1k2(k_1, k_2) = ϕ_1s(k_1) * ϕ_1s(k_2)

    A_kω_Q = single_trion_arpes_signature(
        Q_point, A_SQ_k1k2, 
        m_c, m_v, E_g, E_B, broadening, k_points, ω_points, 
        map(kx -> SVector{2, Float64}([kx, 0.0]), kx_points))     
    
    p = heatmap(kx_points, ω_points, A_kω_Q', 
        c = arpes_colormap(transparency_gradience),
        colorbar = false, 
        legend_foreground_color = :transparent,
        legend = :bottomright,
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
        
    kx_near_valley = LinRange(-0.5, 0.5, 100)
    k_points_near_valley = map(x -> [x, 0.0], kx_near_valley)
    # Position of intensity peak: the brightest spot in one signature, 
    # with varying P_T 
    plot!(p, kx_near_valley, 
        eV * norm.(k_points_near_valley).^2 * (m_T / (2 * m_c^2)) 
        .- eV * norm.(k_points_near_valley).^2 * (m_v / m_c^2)
        .+ E_B .+ E_g,
        label = raw"Varying $\mathbf{P}_\mathrm{T}$, $\mathbf{k}_1 = \mathbf{k}_2 = 0$",
        linestyle = :dot,
        c = colorant"firebrick2")
    # Dispersion relation when P_T is fixed to a constant 
    plot!(p, kx_near_valley, 
        ϵ_v.([(m_c + m_v) / m_T * Q_point] .- k_points_near_valley) 
        .- eV * norm(Q_point)^2 * (m_v / 2m_T^2) .+ eV * norm(Q_point)^2 / 2m_T .+ E_B .+ E_g,
        label = raw"fixed $\mathbf{P}_\mathrm{T}$, $\mathbf{k}_1 = 0$ or $\mathbf{k}_2 = 0$",
        linestyle = :dot, 
        c = colorant"aqua")
    plot!(p, kx_near_valley, 
        2 * ϵ_v.(([Q_point] .- k_points_near_valley) ./ 2) 
        .+ eV * norm(Q_point)^2 / 2m_T .+ E_B .+ E_g,
        label = raw"fixed $\mathbf{P}_\mathrm{T}$, $\mathbf{k}_1 = \mathbf{k}_2 $",
        linestyle = :dot, 
        c = colorant"deepskyblue2")
    
    ylims!(p, (-0.8, 1.8))
    xlims!(p, (-0.8, 0.8))


    cb = heatmap(
        [1 0; 0 1], 
        clims=(0,1), 
        framestyle=:none, 
        c=arpes_colormap(white_gradience), 
        cbar=true, 
        lims=(-1,0)
    )

    p = custom_colorbar(p, cb)
    savefig(p, "trion-P-$(Q_point[1]).png")
end

# Locating the peaks in the trion wave function  
let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_v = 0.8184,
    m_c = 0.4794,
    E_g = 0.7830,
    ϵ = 6.4, 
    E_B = -0.1,
    w = 0.7,
    β = 10,
    broadening(x) = exp(- σ^2 * x^2),
    kx_points = LinRange(-1.2, 1.2, 100),
    ky_points = LinRange(-1.2, 1.2, 100),
    dkx = step(kx_points), 
    dky = step(ky_points),
    dk = dkx * dky,
    k_points = map(t -> SVector{2, Float64}(collect(t)), 
        collect(Iterators.product(kx_points, ky_points))
    ),
    ω_points = LinRange(0, 1.5, 500),
    Q_point  = SVector{2, Float64}([0.5, 0.0]) 

    ϵ_v(k) =  - norm(k)^2 / 2m_v * eV
    ϵ_c(k) =    norm(k)^2 / 2m_c * eV + E_g
    
    m_T = 2m_v + m_c

    μ = m_c * m_v / (m_c + m_v)
    a_ex = 0.529 * ϵ / μ
    ϕ_1s(k) = 1 / (1 + norm(k)^2 * a_ex^2 / 4)^2
    A_SQ_k1k2(k_1, k_2) = ϕ_1s(k_1) * ϕ_1s(k_2)
    
    # Here as a demonstration, we set k_1 to [k_x, 0], and k_2 to [0, k_2]
    heatmap(
        kx_points, ky_points, 
        # The physical meaning of [kx, ky] in k_points
        # now is k_1 = [k_x, 0], k_2 = [0, k_y]
        map(k_points) do (kx, ky)
            abs(A_SQ_k1k2(kx, ky))^2
        end', 
        aspect_ratio = :equal, 
        xlims = (minimum(kx_points) - dkx / 2, maximum(kx_points) + dkx / 2)
    )
end