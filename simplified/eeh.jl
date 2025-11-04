include("model-data.jl")
include("units.jl")
include("qp-bands.jl")
using ThreadTools

"""
Momentum conversion when the electron detected in ARPES comes from peak 1.
Note that in this case k_1 is known from k,
and it's k_2 that needs further specification.
"""
function momentum_calc_eeh_e1(
    trion::Intervalley2DChandraTrion,
    P  :: SVector{2, Float64},
    k  :: SVector{2, Float64},
    k_2:: SVector{2, Float64}, 
)
    m_e = trion.m_e
    m_h = trion.m_h
    M = m_h + 2m_e # This is different from the total mass of the positive trion!
    w = trion.w

    k_e1 = k
    k_1 = k_e1 + m_e / M * w - m_e / M * P
    k_e2 = k_2 + (m_e + m_h) / M * w + m_e / M * P

    k_h = P - k_e1 - k_e2

    (
        k_h = k_h,
        k_e1 = k_e1,
        k_e2 = k_e2,
        P = P,
        k_1 = k_1,
        k_2 = k_2,
    )
end

# Momentum conversion when the electron detected in ARPES comes from peak 2.
# Note that in this case k_2 is known from k,
# and it's k_1 that needs further specification.
function momentum_calc_eeh_e2(
    trion::Intervalley2DChandraTrion,
    P  :: SVector{2, Float64},
    k  :: SVector{2, Float64},
    k_1:: SVector{2, Float64}, 
)
    m_e = trion.m_e
    m_h = trion.m_h
    M = m_h + 2m_e # This is different from the total mass of the positive trion!
    w = trion.w

    k_e2 = k
    k_2 = k_e2 - (m_e + m_h) / M * w - m_e / M * P
    k_e1 = k_1 - m_e / M * w + m_e / M * P

    k_h = P - k_e1 - k_e2

    (
        k_h = k_h,
        k_e1 = k_e1,
        k_e2 = k_e2,
        P = P,
        k_1 = k_1,
        k_2 = k_2,
    )
end

"""
Similar to `momentum_calc_eeh_e1` but calculates all the momentum variables 
from P, ke1, ke2.
"""
function momentum_calc_eeh_rel(
    trion::Intervalley2DChandraTrion,
    P  :: SVector{2, Float64},
    k_e1:: SVector{2, Float64},
    k_e2:: SVector{2, Float64},
)
    m_e = trion.m_e
    m_h = trion.m_h
    M = m_h + 2m_e # This is different from the total mass of the positive trion!
    w = trion.w

    k_1 = k_e1 + m_e / M * (w - P)
    k_2 = k_e2 - (m_e + m_h) / M * w - m_e / M * P 
    k_h = P - k_e1 - k_e2

    (
        k_h = k_h,
        k_e1 = k_e1,
        k_e2 = k_e2,
        P = P,
        k_1 = k_1,
        k_2 = k_2,
    )
end

function E_residue_eeh_e1(
    trion::Intervalley2DChandraTrion,
    k_e2::SVector{2, Float64},
    k_h::SVector{2, Float64}
)
    E_c2(trion, k_e2) + E_v1(trion, k_h)
end

function E_residue_eeh_e2(
    trion::Intervalley2DChandraTrion,
    k_e1::SVector{2, Float64},
    k_h::SVector{2, Float64}
)
    E_c1(trion, k_e1) + E_v1(trion, k_h)
end

function E_trion_eeh(
    trion::Intervalley2DChandraTrion,
    P::SVector{2, Float64}
)
    EB = trion.E_B
    Eg = trion.E_g
    w = trion.w
    m_h = trion.m_h
    m_e = trion.m_e

    M = 2m_e + m_h
    inv_eV * (P - w)' * (P - w) / 2M + 2Eg - EB
end

function trion_ARPES_eeh(
    trion::Intervalley2DChandraTrion,
    P::SVector{2, Float64},
    kgrid::Vector{SVector{2, Float64}},
    kpath::Vector{SVector{2, Float64}},
    freq_list::LinRange{Float64, Int64},
    wfn,
    broaden,
    E_residue_correction::Float64 = 0.0
)
    # ARPES signature from peak 1,
    # where we move k_2 freely
    A_ke1 = tmap(Iterators.product(kpath, freq_list)) do (k, ω)
        sum(kgrid) do k_2
            momentum_set = momentum_calc_eeh_e1(trion, P, k, k_2)
            k_e1 = momentum_set.k_e1
            k_e2 = momentum_set.k_e2
            k_1 = momentum_set.k_1
            k_h = momentum_set.k_h
            broaden(ω - E_trion_eeh(trion, P) + E_residue_eeh_e1(trion, k_e2, k_h) + E_residue_correction) * abs(wfn(k_1, k_2))^2
        end
    end
    
    # ARPES signature from peak 2,
    # where we move k_1 freely.
    A_ke2 = tmap(Iterators.product(kpath, freq_list)) do (k, ω)
        sum(kgrid) do k_1
            momentum_set = momentum_calc_eeh_e2(trion, P, k, k_1)
            k_e1 = momentum_set.k_e1
            k_e2 = momentum_set.k_e2
            k_2 = momentum_set.k_2
            k_h = momentum_set.k_h
            broaden(ω - E_trion_eeh(trion, P) + E_residue_eeh_e2(trion, k_e1, k_h) + E_residue_correction) * abs(wfn(k_1, k_2))^2
        end
    end
    
    A_ke1 + A_ke2
end
