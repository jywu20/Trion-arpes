using ProgressMeter
include("exciton-solver.jl")

m_h = 0.21
m_e = 0.37
E_g = 2.84
ϵ = 6.4 
E_B = -0.1
hexagonal_edge = 4π / (3 * 3.144817974)
w = hexagonal_edge
β = 1
a = 10.3
b = 25.2
kx_points = LinRange(-0.35, 1.7, 250)
ky_points = LinRange(-0.35, 1.7, 250)
dkx = step(kx_points)
dky = step(ky_points)
dk = dkx * dky
k_points = map(t -> SVector{2, Float64}(collect(t)), 
    collect(Iterators.product(kx_points, ky_points))
)
ω_points = LinRange(-1, 3, 500)
P_x = hexagonal_edge * 0.8
P_T  = SVector{2, Float64}([P_x, 0.0]) 
k_e = SVector{2, Float64}([-m_e / M * (w - P_x) + 0.1, 0.0])

ϵ_v2(k) =  - norm(k - w)^2 / 2m_h * inv_eV
ϵ_v1(k) =  - norm(k)^2     / 2m_h * inv_eV
ϵ_c1(k) =    norm(k)^2     / 2m_e * inv_eV + E_g
ϵ_c2(k) =    norm(k - w)^2 / 2m_e * inv_eV + E_g
E_SP(P) = inv_eV * norm(P .- w)^2 / 2M .+ E_B .+ 2E_g
