# Illustration of the profiles of broadening functions

using CairoMakie
include("broadening.jl")

let f = Figure()
    Es = LinRange(-5, 5, 100)
    ax = Axis(f[1, 1])

    for σ in [2.0, 10.0, 20.0]
        broaden = gaussian_broadening(σ)
        lines!(ax, Es, broaden.(Es), label="$σ")
    end
    axislegend(ax)
        
    f
end
