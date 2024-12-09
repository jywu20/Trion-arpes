include("model-data.jl")
include("units.jl")
using ThreadTools

"""
Momentum transform for the positive trion mode, containing two holes and one electron.
"""
function rel2abs_ehh(
    trion::Intervalley2DChandraTrion, 
    k_1::SVector{2, Float64}, 
    k_2::SVector{2, Float64},
    P  ::SVector{2, Float64}
)
    m_e = trion.m_e
    m_h = trion.m_h
    M = 2m_h + m_e
    w = trion.w

    k_e  = (P - w) * m_e/M - k_1 - k_2
    k_h1 = k_1 - m_h/M * w + m_h/M * P
    k_h2 = k_2 + (m_e+m_h)/M * w + m_h/M * P

    k_e, k_h1, k_h2
end

function momentum_calc_ehh(
    trion::Intervalley2DChandraTrion,
    P  :: SVector{2, Float64},
    k  :: SVector{2, Float64},
    k_1:: SVector{2, Float64}, 
)
    m_e = trion.m_e
    m_h = trion.m_h
    M = 2m_h + m_e
    w = trion.w

    k_e = k
    k_2 = (P - w) * m_e/M - k_1 - k_e

    k_h1 = k_1 - m_h/M * w + m_h/M * P
    k_h2 = k_2 + (m_e+m_h)/M * w + m_h/M * P

    (
        k_e = k_e,
        k_h1 = k_h1,
        k_h2 = k_h2,
        P = P,
        k_1 = k_1,
        k_2 = k_2,
    )
end

function E_residue_ehh(
    trion::Intervalley2DChandraTrion,
    k_h1::SVector{2, Float64},
    k_h2::SVector{2, Float64}
)
    E_v1(trion, k_h1) + E_v2(trion, k_h2)
end

function E_trion(
    trion::Intervalley2DChandraTrion,
    P::SVector{2, Float64}
)
    EB = trion.E_B
    Eg = trion.E_g
    w = trion.w
    m_h = trion.m_h
    m_e = trion.m_e

    M = 2m_h + m_e
    inv_eV * (P - w)' * (P - w) / 2M + Eg - EB
end

function trion_ARPES_ehh(
    trion::Intervalley2DChandraTrion,
    P::SVector{2, Float64},
    kgrid::Vector{SVector{2, Float64}},
    kpath::Vector{SVector{2, Float64}},
    freq_list::LinRange{Float64, Int64}
)
    
end