include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("ehh.jl")
include("wfn.jl")
include("arpes.jl")
using LinearAlgebra

m_h = 0.21
m_e = 0.37
E_g = 2.84
w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])
trion = Intervalley2DChandraTrion(
    m_h = m_h,
    m_e = m_e,
    w = w,
    E_g = E_g,
    E_B = 0.75, # Binding energy for the ehh trion mode
    a = 10.3,
    b = 25.2
)
M = 2m_h + m_e
E_B = 0.75

# The trion momentum is set to be w
P_ratio = 1.2
P = P_ratio * w

kx_list = LinRange(-0.5, 0.5, 200)
k1_list = [SA[kx, 0.0] for kx in kx_list]
k1_grid = reshape(map(Iterators.product(kx_list, kx_list)) do (k_x, k_y)
    SA[k_x, k_y]
end, length(k1_list)^2)

k = SA[0.2, -0.3]

dispersion_k_zero = map(k1_list) do k
    k_1 = SA[0.0, 0.0]
    momentum_set = momentum_calc_ehh(trion, P, k, k_1)
    k_h1 = momentum_set.k_h1
    k_h2 = momentum_set.k_h2 
    k_2  = momentum_set.k_2
    E_trion_ehh(trion, P) - E_residue_ehh(trion, k_h1, k_h2)
end

dispersion_k_zero_analytic = map(k1_list) do k
    inv_eV * (
        - norm(k - (m_e + m_h)/M*(P-w))^2 / 2m_h  
        + norm(P-w)^2 * (m_e + m_h) / 2M^2
    ) + E_g - E_B
end

println(norm(dispersion_k_zero - dispersion_k_zero_analytic))

dispersion_k_equal = map(k1_list) do k
    M = 2trion.m_h + trion.m_e
    k_1 = ((P - w) * trion.m_e / M - k) / 2
    k_2 = k_1
    momentum_set = momentum_calc_ehh(trion, P, k, k_1)
    k_h1 = momentum_set.k_h1
    k_h2 = momentum_set.k_h2 
    E_trion_ehh(trion, P) - E_residue_ehh(trion, k_h1, k_h2)
end

dispersion_k_equal_analytic = map(k1_list) do k
    inv_eV * (-norm(P - w - k)^2/4m_h + norm(P-w)^2/2M) + E_g - E_B
end

println(norm(dispersion_k_equal - dispersion_k_equal_analytic))