# Plot the ARPES heatmap;

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("exciton.jl")
include("ehh.jl")
include("wfn.jl")
include("arpes.jl")
include("overlap.jl")
include("exciton-band-gamma.jl")
include("exciton-band-k.jl")
include("bgw-crystal.jl")
using CairoMakie
using LaTeXStrings
using Colors

k_K = SVector{3, Float64}(0.333333333333333, 0.3333333333333, 0)
rk, Avck = read_ex_wfc("../../MoS2/MoS2/4-absorption-120-no-sym-gw/eigenvectors.h5", k_K)
k_K_real = read_B("../../MoS2/MoS2/4-absorption-120-no-sym-gw/eigenvectors.h5") * k_K

function fetch_S(Avck::Array{ComplexF64, 4}, iS::Int)
    Avck_S = Avck[:, :, :, iS]
    Avck_S = Avck_S[1, 1, :] + Avck_S[1, 2, :]
    Avck_S
end

A1s_like = fetch_S(Avck, 2)
A2p_like_1 = fetch_S(Avck, 8)
A2p_like_2 = fetch_S(Avck, 6)
A2s_like = fetch_S(Avck, 10)

phase_factor = A1s_like ./ abs.(A1s_like)

let f = Figure(size=(900, 400))
    Avck_plot = A2p_like_1 ./ phase_factor
    Avck_plot ./= Avck_plot[1] / abs.(Avck_plot[1])
    intensity_limit = maximum(abs.(Avck_plot))
    
    Avck_plot_Re = real.(Avck_plot)
    ax = Axis(f[1, 1], title="$(norm(Avck_plot_Re))")
    scatter!(ax, rk[1, :], rk[2, :], color=Avck_plot_Re, 
        colormap=cgrad(:balance),
        colorrange=(-intensity_limit, intensity_limit),
    )
    colsize!(f.layout, 1, Aspect(1, 1))
    hidedecorations!(ax, ticklabels = false, ticks = false)
    
    Avck_plot_Im = imag.(Avck_plot)
    ax = Axis(f[1, 2], title="$(norm(Avck_plot_Im))")
    scatter!(ax, rk[1, :], rk[2, :], color=Avck_plot_Im, 
        colormap=cgrad(:balance),
        colorrange=(-intensity_limit, intensity_limit),
    )
    colsize!(f.layout, 2, Aspect(1, 1))
    hidedecorations!(ax, ticklabels = false, ticks = false)

    f
end
