using CairoMakie
include("wfn.jl")
include("eeh.jl")
include("qp-bands.jl")
include("broadening.jl")
include("units.jl")

kx_list = LinRange(-0.5, 0.5, 200)
ky_list = LinRange(-0.5, 0.5, 200)
ω_list = LinRange(0, 3, 100)

k2_list = Iterators.product(kx_list, ky_list)
w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0]) 

trion = Intervalley2DChandraTrion(
    m_h = 0.21,
    m_e = 0.37,
    w = w,
    E_g = 2.84,
    E_B = 0.76,
    a = 10.3,
    b = 25.2
)
A_k1k2 = wfn(trion)

k = SVector{2, Float64}([-0.2, 0.0])
P = 0.8w + SA[0.0, 0.1]

# Peak 1 is being detected, and we change k_2
Asq_k2 = map(k2_list) do (k_x, k_y)
    k_2 = SVector{2, Float64}([k_x, k_y])
    k_1 = momentum_calc_eeh_e1(trion, P, k, k_2).k_1
    abs2(A_k1k2(k_1, k_2))
end

f = Figure(size=(400, 800)) 
aspect_ratio = (maximum(kx_list) - minimum(kx_list)) / (maximum(ky_list) - minimum(ky_list))

ax = Axis(f[1, 1], aspect=aspect_ratio, title="Wave function structural factor")
heatmap!(ax, kx_list, ky_list, Asq_k2)

let mom_report = momentum_calc_eeh_e1(trion, P, k, SA[0.0, 0.0])
    k_2 = mom_report.k_2
    vlines!(ax, [k_2[1]], color=:white, linewidth=0.5)
    hlines!(ax, [k_2[2]], color=:white, linewidth=0.5)
end

ω = 2.0
σ = 20.0fs

broaden = gaussian_broadening(σ)

ax = Axis(f[2, 1], aspect=aspect_ratio, title="Dispersion relation constraint")
E_residue_k2 = map(k2_list) do (k_x, k_y)
    k_2 = SVector{2, Float64}([k_x, k_y])
    k_e2 = momentum_calc_eeh_e1(trion, P, k, k_2).k_e2
    k_h = momentum_calc_eeh_e1(trion, P, k, k_2).k_h
    broaden(ω - E_trion_eeh(trion, P) + E_residue_eeh_e1(trion, k_e2, k_h)) 
end

heatmap!(ax, kx_list, ky_list, E_residue_k2, 
    colorrange=(minimum(E_residue_k2), maximum(E_residue_k2)))

let mom_report = momentum_calc_eeh_e1(trion, P, k, SA[0.0, 0.0])
    k_2 = mom_report.k_2
    vlines!(ax, [k_2[1]], color=:white, linewidth=0.5)
    hlines!(ax, [k_2[2]], color=:white, linewidth=0.5)
end

Colorbar(f[2, 2],  colorrange=(minimum(E_residue_k2), maximum(E_residue_k2)))

ax = Axis(f[3, 1], aspect=aspect_ratio, title="Contribution from one k_2")
Asq_k2_with_constraints = E_residue_k2 .* Asq_k2
heatmap!(ax, kx_list, ky_list, Asq_k2_with_constraints)

let mom_report = momentum_calc_eeh_e1(trion, P, k, SA[0.0, 0.0])
    k_2 = mom_report.k_2
    vlines!(ax, [k_2[1]], color=:white, linewidth=0.5)
    hlines!(ax, [k_2[2]], color=:white, linewidth=0.5)
end

Colorbar(f[3, 2],  colorrange=(minimum(Asq_k2_with_constraints), maximum(Asq_k2_with_constraints)))

save("eeh-ring-and-peak-contribution-from-k2-example.png", f)

f