include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("exciton.jl")
include("ehh.jl")
include("wfn.jl")
include("arpes.jl")
include("overlap.jl")
using CairoMakie
using LaTeXStrings
using Colors

rk, Avck = read_ex_wfc("../../MoS2/MoS2/4-absorption-120-no-sym-gw/eigenvectors.h5", SVector{3, Float64}(0.333333333333333, 0.3333333333333, 0))

Avck_1s = Avck[1, 1, :, 2] + Avck[1, 2, :, 2]
Avck_1s = abs.(Avck_1s)

let f = Figure(size=(800, 400))
    ax = Axis(f[1, 1])
    scatter!(ax, rk[1, :], rk[2, :], color=Avck_1s, 
        colormap=reverse(cgrad(:grayC)),
        colorrange=(0, 1),
    )
    colsize!(f.layout, 1, Aspect(1, 1))

    ax = Axis(f[1, 2])
    rk_lengths = map(1 : size(rk)[2]) do ik
        norm(rk[:, ik])
    end
    
    scatter!(ax, rk_lengths, Avck_1s)
    
    a = 10.4
    rk_lengths_sample = LinRange(0, maximum(rk_lengths), 400)
    lines!(ax, rk_lengths_sample, map(rk_lengths_sample) do k
        1 / (1 + k^2 * a^2)^1.5 * maximum(Avck_1s)
    end)

    colsize!(f.layout, 2, Aspect(1, 1))

    save("exciton-radius.png", f)
    f
end