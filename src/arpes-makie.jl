using CairoMakie
using MakieTeX
using LinearAlgebra
using Colors
using ColorSchemes
using LaTeXStrings

function arpes_colormap(gradient)
    rainbow = cgrad(:rainbow)
    starting_region = gradient(first(rainbow), steps = 50) 
    middle_region = rainbow[LinRange(0, 1, 500)]
    end_region = reverse(gradient(last(rainbow), steps = 400))
    
    vcat(starting_region, middle_region, end_region)
end

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

"""
The main purpose of this function is to avoid the broadening function 
depending on a global variable σ,
whose type is variable, which makes it slow.
"""
function gaussian_broadening(σ::Float64)
    broadening(x::Float64) = @fastmath exp(- σ^2 * x^2)
    broadening
end