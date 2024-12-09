# A band plot, showing that the broadening function 
# and the definition of the dispersion relations is not wrong;
# possibly as a benchmark with the first-principles MoS2 bands

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
using CairoMakie

w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])
trion = Intervalley2DChandraTrion(
    m_h = 0.54,
    m_e = 0.60,
    w = w,
    E_g = 2.84,
    E_B = 0.1,
    a = 10.3,
    b = 25.2
)

E_c_GW = Float64[]
E_v_GW = Float64[]

iK, iK′ = open("MoS2-K-Kp.dat") do f
    iK = 0
    iK′ = 0
    readline(f)
    readline(f)
    while !eof(f)
        tokens = split(readline(f))
        ik = parse(Int, tokens[1])
        kx = parse(Float64, tokens[2])
        push!(E_v_GW, parse(Float64, tokens[3]))
        push!(E_c_GW, parse(Float64, tokens[4]))
        if abs(kx - 0.33333) < 1e-4
            iK = ik
        elseif abs(kx - 0.66666) < 1e-4
            iK′ = ik
        end
    end
    iK, iK′
end
k_list_GW = eachindex(E_c_GW)
k_list_GW = k_list_GW .- iK
k_list_GW = k_list_GW .* (w_side / (iK′ - iK))
E_v_GW .-= maximum(E_v_GW)
E_c_GW .-= minimum(E_c_GW) .- trion.E_g

f = Figure()

kx_list = LinRange(-0.35, 0.35, 250)
k_list = [SA[kx, 0.0] for kx in kx_list]
ω_list = LinRange(-2.2, 4.2, 100)

E_c1_list = E_c1.([trion], k_list)
E_v1_list = E_v1.([trion], k_list)

broaden = gaussian_broadening(10.0)

ax = Axis(f[1, 1])
heatmap!(ax, kx_list, ω_list, map(Iterators.product(kx_list, ω_list)) do (k_x, ω)
    k_e = SA[k_x, 0.0]
    broaden(ω - E_c1(trion, k_e))
end)
lines!(ax, kx_list, E_c1_list)
scatter!(ax, k_list_GW, E_c_GW, color=:green)
xlims!(ax, (minimum(kx_list), maximum(kx_list)))
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

ax = Axis(f[1, 2])
heatmap!(ax, kx_list, ω_list, map(Iterators.product(kx_list, ω_list)) do (k_x, ω)
    k_e = SA[k_x, 0.0]
    broaden(ω + E_v1(trion, k_e))
end)
lines!(ax, kx_list, -E_v1_list)
scatter!(ax, k_list_GW, E_v_GW, color=:green)
xlims!(ax, (minimum(kx_list), maximum(kx_list)))
hidedecorations!(ax, ticks = false, ticklabels = false, label = false)

save("effective-mass-benchmark.png", f)