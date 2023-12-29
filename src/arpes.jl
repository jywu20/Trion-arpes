using Plots
using LinearAlgebra
using Colors
using ColorSchemes
using LaTeXStrings

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

