using StaticArrays
using ThreadTools

include("real-space.jl")
include("two-band.jl")


function trionarpes_e2h1h1_def(
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

function trionarpes_e2h1h1_thread(
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

function trionarpes_e1h1h2_thread(
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
            ϵ_v_kh2 = - inv_eV * (k_h2 - w)' * (k_h2 - w) / 2m_h

            δ_factor = broadening(ω - E_S_PT - ϵ_v_kh1 - ϵ_v_kh2 - E_B - E_g)
            
            k_1 = k_h1 + m_h / M * w - m_h / M * P_T
            k_2 = k_h2 - (m_h + m_e) / M * w - m_h / M * P_T
            trion_structure_factor = abs(A_SQ_k1k2(k_1, k_2))^2
            
            δ_factor * trion_structure_factor
        end
    end
    
    A_kω
end

function trionarpes_e1e2h1_thread(
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
    
    A_kω_1 = tmap(heatmap_points) do (k_e1, ω)
        sum(k_points) do k_e2
            k_h = P_T - k_e1 - k_e2
            ϵ_c2 =   inv_eV * (k_e2 - w)' * (k_e2 - w) / 2m_e
            ϵ_v  = - inv_eV * k_h' * k_h / 2m_h

            δ_factor = broadening(ω - E_S_PT - ϵ_c2 - ϵ_v - E_B - E_g)
            
            k_1 = k_e1 + m_e / M * w - m_e / M * P_T
            k_2 = k_e2 - (m_e + m_h) / M * w - m_e / M * P_T
            trion_structure_factor = abs(A_SQ_k1k2(k_1, k_2))^2
            
            δ_factor * trion_structure_factor
        end
    end

    A_kω_2 = tmap(heatmap_points) do (k_e2, ω)
        sum(k_points) do k_e1
            k_h = P_T - k_e1 - k_e2
            ϵ_c1 =   inv_eV * k_e1' * k_e1 / 2m_e
            ϵ_v  = - inv_eV * k_h' * k_h / 2m_h

            δ_factor = broadening(ω - E_S_PT - ϵ_c1 - ϵ_v - E_B - E_g)
            
            k_1 = k_e1 + m_e / M * w - m_e / M * P_T
            k_2 = k_e2 - (m_e + m_h) / M * w - m_e / M * P_T
            trion_structure_factor = abs(A_SQ_k1k2(k_1, k_2))^2
            
            δ_factor * trion_structure_factor
        end
    end
    
    A_kω_1 + A_kω_2
end

struct IndirectTwoBandMat2D 
    ham::IndirectTwoBandModel2D
    dielectric::Dielectric2D
    # Possible new fields: Berry curvature, dipole, etc.
end

function phi1s_def(trion_spec::IndirectTwoBandMat2D)
    ham = trion_spec.ham
    dielectric = trion_spec.dielectric

    m_c = ham.m_e
    m_v = ham.m_h
    ϵ = dielectric.ϵ

    μ = m_c * m_v / (m_c + m_v)
    a_ex = a_Bohr * ϵ / μ

    A_SQ_k1k2(k_1, k_2) = ϕ_1s_def(k_1, a_ex) * ϕ_1s_def(k_2, a_ex)
    
    A_SQ_k1k2
end

function phi1s(trion_spec::IndirectTwoBandMat2D)
    ham = trion_spec.ham
    dielectric = trion_spec.dielectric

    m_c = ham.m_e
    m_v = ham.m_h
    ϵ = dielectric.ϵ

    μ = m_c * m_v / (m_c + m_v)
    a_ex = a_Bohr * ϵ / μ

    @inline A_SQ_k1k2(k_1, k_2) = ϕ_1s(k_1, a_ex) * ϕ_1s(k_2, a_ex)
    
    A_SQ_k1k2
end

function phi1sa1sb(trion_spec::IndirectTwoBandMat2D, a::Float64, b::Float64)
    ham = trion_spec.ham
    dielectric = trion_spec.dielectric

    m_c = ham.m_e
    m_v = ham.m_h
    ϵ = dielectric.ϵ

    μ = m_c * m_v / (m_c + m_v)
    a_ex = a_Bohr * ϵ / μ

    @inline ϕ_1s_a(k) = ϕ_1s(k, a)  
    @inline ϕ_1s_b(k) = ϕ_1s(k, b)
    @inline A_SQ_k1k2(k_1, k_2) = 1/sqrt(2) * (
        ϕ_1s_a(k_1) * ϕ_1s_b(k_2) + ϕ_1s_a(k_2) * ϕ_1s_b(k_1)
    )
    
    A_SQ_k1k2
end
