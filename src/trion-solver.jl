using StaticArrays
using ThreadTools

const fs = 1.5193    
const inv_eV = 3.80998

struct IndirectTwoBandModel2D
    m_e :: Float64     # Conduction band effective mass
    m_h :: Float64     # Valance band effective mass (positive)
    E_g :: Float64     # Single particle band gap
    w   :: SVector{2, Float64} # k_{valley} - k_{peak}
end

struct Dielectric2D 
    ϵ :: Float64
    # TODO
end


function single_trion_arpes_signature_def(
    ham::IndirectTwoBandModel2D, 
    E_B,                                 # Binding energy in eV
    A_SQ_k1k2,                           # The internal part of the trion wave function  
                                         # The trion band index S and momentum Q should be pre-determined
    P_T::SVector{2, Float64},            # Total momentum of the trion  
    k_points,       # The 1BZ sampling 
    ω_points,       # The frequency sampling, in eV 
    k_points_ARPES, # k points included in the ARPES plotting
    broadening,     # Broadening function
)

    m_e = ham.m_e
    m_h = ham.m_h
    E_g = ham.E_g
    M = 2m_h + m_e
    w = ham.w

    E_S_PT = inv_eV * norm(P_T - w)^2 / 2M

    heatmap_points = Iterators.product(k_points_ARPES, ω_points) 

    A_kω = map(heatmap_points) do (k_e, ω)
        sum(k_points) do k_h1
            k_h2 = P_T - k_e - k_h1
            ϵ_v_kh1 = - inv_eV * norm(k_h1)^2 / 2m_h
            ϵ_v_kh2 = - inv_eV * norm(k_h2)^2 / 2m_h

            δ_factor = broadening(ω - E_S_PT - ϵ_v_kh1 - ϵ_v_kh2 - E_B - E_g)
            
            k_1 = k_h1 + m_h / M * w - m_h / M * P_T
            k_2 = k_h2 + m_h / M * w - m_h / M * P_T
            trion_structure_factor = abs(A_SQ_k1k2(k_1, k_2))^2
            
            δ_factor * trion_structure_factor
        end
    end
    
    A_kω
end

single_trion_arpes_signature = single_trion_arpes_signature_def

function single_trion_arpes_signature_thread(
    ham::IndirectTwoBandModel2D, 
    E_B,                                 # Binding energy in eV
    A_SQ_k1k2,                           # The internal part of the trion wave function  
    P_T::SVector{2, Float64},            # Total momentum of the trion  
    k_points,       # The 1BZ sampling 
    ω_points,       # The frequency sampling, in eV 
    k_points_ARPES, # k points included in the ARPES plotting
    broadening,     # Broadening function
)

    m_e = ham.m_e
    m_h = ham.m_h
    E_g = ham.E_g
    M = 2m_h + m_e
    w = ham.w

    E_S_PT = inv_eV * (P_T - w)' * (P_T - w) / 2M

    heatmap_points = Iterators.product(k_points_ARPES, ω_points) 
    
    A_kω = tmap(heatmap_points) do (k_e, ω)
        sum(k_points) do k_h1
            k_h2 = P_T - k_e - k_h1
            ϵ_v_kh1 = - inv_eV * k_h1' * k_h1 / 2m_h
            ϵ_v_kh2 = - inv_eV * k_h2' * k_h2 / 2m_h

            δ_factor = broadening(ω - E_S_PT - ϵ_v_kh1 - ϵ_v_kh2 - E_B - E_g)
            
            k_1 = k_h1 + m_h / M * w - m_h / M * P_T
            k_2 = k_h2 + m_h / M * w - m_h / M * P_T
            trion_structure_factor = abs(A_SQ_k1k2(k_1, k_2))^2
            
            δ_factor * trion_structure_factor
        end
    end
    
    A_kω
end

ϕ_1s_def(k::SVector{2, Float64}, a_ex::Float64) = 1 / (1 + norm(k)^2 * a_ex^2 / 4)^2
@inline ϕ_1s(k::SVector{2, Float64}, a_ex::Float64) = @fastmath 1 / (1 + k' * k * a_ex^2 / 4)^2

function ground_state_1s_def(ham::IndirectTwoBandModel2D, dielectric::Dielectric2D)
    m_c = ham.m_e
    m_v = ham.m_h
    ϵ = dielectric.ϵ

    μ = m_c * m_v / (m_c + m_v)
    a_ex = 0.529 * ϵ / μ

    A_SQ_k1k2(k_1, k_2) = ϕ_1s_def(k_1, a_ex) * ϕ_1s_def(k_2, a_ex)
    
    A_SQ_k1k2
end

function ground_state_1s(ham::IndirectTwoBandModel2D, dielectric::Dielectric2D)
    m_c = ham.m_e
    m_v = ham.m_h
    ϵ = dielectric.ϵ

    μ = m_c * m_v / (m_c + m_v)
    a_ex = 0.529 * ϵ / μ

    @inline A_SQ_k1k2(k_1, k_2) = ϕ_1s(k_1, a_ex) * ϕ_1s(k_2, a_ex)
    
    A_SQ_k1k2
end