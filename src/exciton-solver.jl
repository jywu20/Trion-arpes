using StaticArrays
using ThreadTools

include("two-band.jl")

function single_exciton_arpes_signature_def(
    ham::IndirectTwoBandModel2D,
    E_B, 

    A_SQ_k,
    Q::SVector{2, Float64},
    ω_points,
    k_points_ARPES,
    broadening,
)
    m_e = ham.m_e
    m_h = ham.m_h
    E_g = ham.E_g
    w = ham.w
    
    M = m_e + m_h
    E_SQ = inv_eV * norm(Q - w)^2 / 2M
    
    heatmap_points = Iterators.product(k_points_ARPES, ω_points)
    A_kω = map(heatmap_points) do (k_e, ω)
        k_h = Q - k_e

        ϵ_v = - inv_eV * norm(k_h)^2 / 2m_h

        δ_factor = broadening(ω - E_SQ - ϵ_v - E_B - E_g)
        
        k = k_e - m_e * Q / M - m_h * w / M
        exciton_structure_factor = abs(A_SQ_k(k))^2
        
        δ_factor * exciton_structure_factor
    end
    
    A_kω
end

struct IndirectTwoBandExciton2D 
    ham::IndirectTwoBandModel2D
    dielectric::Dielectric2D
end

function ground_state_1s(exciton_spec::IndirectTwoBandExciton2D)
    ham = exciton_spec.ham
    dielectric = exciton_spec.dielectric

    m_c = ham.m_e
    m_v = ham.m_h
    ϵ = dielectric.ϵ

    μ = m_c * m_v / (m_c + m_v)
    a_ex = 0.529 * ϵ / μ

    x -> ϕ_1s(x, a_ex)    
end