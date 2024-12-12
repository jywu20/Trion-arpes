using CairoMakie
include("wfn.jl")
include("ehh.jl")
include("qp-bands.jl")
include("broadening.jl")
include("units.jl")

kx_list = LinRange(-0.5, 0.5, 200)
ky_list = LinRange(-0.5, 0.5, 200)
ω_list = LinRange(0, 3, 100)
#kx_list = LinRange(-1, 1, 100)
#ky_list = LinRange(-1, 1, 100)

k1_list = Iterators.product(kx_list, ky_list)
w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0]) + SA[0, 0.1]

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

k = SVector{2, Float64}([0.3, -0.1])
P = 1.2w


Asq = map(k1_list) do (k_x, k_y)
    k_1 = SVector{2, Float64}([k_x, k_y])
    k_2 = momentum_calc_ehh(trion, P, k, k_1).k_2
    abs2(A_k1k2(k_1, k_2))
end

f = Figure(size=(400, 800)) 
aspect_ratio = (maximum(kx_list) - minimum(kx_list)) / (maximum(ky_list) - minimum(ky_list))

ax = Axis(f[1, 1], aspect=aspect_ratio, title="Wave function structural factor")
heatmap!(ax, kx_list, ky_list, Asq)

let mom_report = momentum_calc_ehh(trion, P, k, SA[0.0, 0.0])
    k_1 = mom_report.k_1
    k_2 = mom_report.k_2
    vlines!(ax, [k_1[1]], color=:white, linewidth=0.5)
    vlines!(ax, [k_2[1]], color=:white, linewidth=0.5)
    hlines!(ax, [k_1[2]], color=:white, linewidth=0.5)
    hlines!(ax, [k_2[2]], color=:white, linewidth=0.5)
end

ω = 2.3
σ = 20.0fs

broaden = gaussian_broadening(σ)

ax = Axis(f[2, 1], aspect=aspect_ratio, title="Dispersion relation constraint")
E_residue = map(k1_list) do (k_x, k_y)
    k_1 = SVector{2, Float64}([k_x, k_y])
    k_h1 = momentum_calc_ehh(trion, P, k, k_1).k_h1
    k_h2 = momentum_calc_ehh(trion, P, k, k_1).k_h2 
    broaden(ω - E_trion_ehh(trion, P) + E_residue_ehh(trion, k_h1, k_h2))
end

heatmap!(ax, kx_list, ky_list, E_residue, 
    colorrange=(minimum(E_residue), maximum(E_residue)))

let mom_report = momentum_calc_ehh(trion, P, k, SA[0.0, 0.0])
    k_1 = mom_report.k_1
    k_2 = mom_report.k_2
    vlines!(ax, [k_1[1]], color=:white, linewidth=0.5)
    vlines!(ax, [k_2[1]], color=:white, linewidth=0.5)
    hlines!(ax, [k_1[2]], color=:white, linewidth=0.5)
    hlines!(ax, [k_2[2]], color=:white, linewidth=0.5)
end

Colorbar(f[2, 2],  colorrange=(minimum(E_residue), maximum(E_residue)))

ax = Axis(f[3, 1], aspect=aspect_ratio, title="Final contribution from one k_1")
Asq_with_constraint = E_residue .* Asq
heatmap!(ax, kx_list, ky_list, Asq_with_constraint)
Colorbar(f[3, 2], 
    colorrange=(minimum(Asq_with_constraint), maximum(Asq_with_constraint)))

let mom_report = momentum_calc_ehh(trion, P, k, SA[0.0, 0.0])
    k_1 = mom_report.k_1
    k_2 = mom_report.k_2
    vlines!(ax, [k_1[1]], color=:white, linewidth=0.5)
    vlines!(ax, [k_2[1]], color=:white, linewidth=0.5)
    hlines!(ax, [k_1[2]], color=:white, linewidth=0.5)
    hlines!(ax, [k_2[2]], color=:white, linewidth=0.5)
end

println(sum(Asq_with_constraint))

save("ehh-ring-and-peak-contribution-example.png", f)

f
