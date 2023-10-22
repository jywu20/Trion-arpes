using Plots
using LinearAlgebra
using ProgressMeter
using Colors
using ColorSchemes

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
    starting_region = gradient(first(rainbow), steps = 100) 
    middle_region = rainbow[LinRange(0, 1, 500)]
    end_region = reverse(gradient(last(rainbow), steps = 400))
    
    PlotUtils.ContinuousColorGradient([starting_region..., middle_region..., end_region...])
end

function custom_colorbar(p, cb)
    l = @layout [a{0.95w} b]
    plot(p, cb, layout = l)
end


# The exciton band index S and momentum Q should be pre-determined
function single_exciton_arpes_signature(ϵ_v, Q, A_SQ, E_SQ, broadening, k_points, ω_points)
    heatmap_points = Iterators.product(k_points, ω_points) 

    δ_factor = map(heatmap_points) do (k, ω)
        broadening(ω - E_SQ - ϵ_v(k - Q))
    end
    
    exciton_structure_factor = map(heatmap_points) do (k, ω)
        abs(A_SQ(k))^2
    end
    
    A_kω = exciton_structure_factor .* δ_factor
    A_kω
end

# Example of ARPES plot 

const fs = 1.5193
eV = 3.80998

p = let σ = 20fs, # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
    m_v = 0.8184,
    m_c = 0.4794,
    E_g = 0.7830,
    ϵ = 6.4, 
    E_B = 0.1,
    w = 0.7,
    Q = w + 0.2,
    β = 100,
    broadening(x) = exp(- σ^2 * x^2),
    k_points = LinRange(-0.5, 1.2, 500),
    ω_points = LinRange(-0.4, 1.2, 500),

    ϵ_v(k) = - norm(k)^2 / 2m_v * eV
    ϵ_c(k) = E_g + norm(k - w)^2 / 2m_c * eV

    M = m_c + m_v
    α = m_c / M
    μ = m_c * m_v / (m_c + m_v)
    a_ex = 0.529 * ϵ / μ
    A_SQ(k) = 1 / (1 + norm(k - (w + α * (Q - w)))^2 * a_ex^2 / 4)^2
    E_SQ = E_g - E_B


    A_kω = single_exciton_arpes_signature(ϵ_v, Q, A_SQ, E_SQ, broadening, k_points, ω_points)
    
    p = heatmap(k_points, ω_points, A_kω', 
        c = arpes_colormap(transparency_gradience),
        colorbar = false) 
    plot!(p, k_points, ϵ_v.(k_points), 
        legend = false, 
        linestyle = :dash)
    plot!(p, k_points, ϵ_c.(k_points), 
        legend = false, 
        linestyle = :dash)
    
    xlims!(p, (-0.5, 1.2))
    ylims!(p, (-0.4, 1.2)) 
    
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

