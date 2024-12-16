include("model-data.jl")
include("wfn.jl")

"""
Note that here `k` is the momentum of the photoelectron,
and not the internal momentum between the electron and the hole.
"""
function momentum_calc(
    exciton::IntraValley2DExciton,
    Q::SVector{2, Float64},
    k::SVector{2, Float64},
)
    m_e = exciton.m_e
    m_h = exciton.m_h
    M = m_e + m_h

    k_inner = k - m_e * Q / M 
    k_e = k
    k_h = Q - k

    (
        k_e = k_e,
        k_h = k_h,
        Q = Q,
        k_i = k_inner, 
    )
end

"""
Note that here `k` is the momentum of the photoelectron,
and not the internal momentum between the electron and the hole.
"""
function momentum_calc(
    exciton::InterValley2DExciton,
    Q::SVector{2, Float64},
    k::SVector{2, Float64},
)
    m_e = exciton.m_e
    m_h = exciton.m_h
    M = m_e + m_h
    w = exciton.w

    k_inner = k - m_e * Q / M - m_h * w / M
    k_e = k
    k_h = Q - k

    (
        k_e = k_e,
        k_h = k_h,
        Q = Q,
        k_i = k_inner, 
    )
end

function E_exciton(
    exciton::IntraValley2DExciton,
    Q::SVector{2, Float64},
)
    M = exciton.m_e + exciton.m_h
    E_B = exciton.E_B
    E_g = exciton.E_g
    inv_eV * Q' * Q / 2M + E_g - E_B
end

function E_exciton(
    exciton::InterValley2DExciton,
    Q::SVector{2, Float64},
)
    M = exciton.m_e + exciton.m_h
    E_B = exciton.E_B
    E_g = exciton.E_g
    w = exciton.w
    inv_eV * (Q - w)' * (Q - w) / 2M + E_g - E_B
end

function exciton_ARPES_eh(
    exciton::IntraValley2DExciton,
    Q::SVector{2, Float64},
    kpath::Vector{SVector{2, Float64}},
    freq_list::LinRange{Float64, Int64},
    wfn,
    broaden
)
    map(Iterators.product(kpath, freq_list)) do (k, ω)
        momentum_set = momentum_calc(exciton, Q, k)
        k_h = momentum_set.k_h
        k_e = momentum_set.k_e
        k_i = momentum_set.k_i
        broaden(ω - E_exciton(exciton, Q) + E_v(exciton, k_h)) * abs2(wfn(k_i))
    end
end

function exciton_ARPES_eh(
    exciton::InterValley2DExciton,
    Q::SVector{2, Float64},
    kpath::Vector{SVector{2, Float64}},
    freq_list::LinRange{Float64, Int64},
    wfn,
    broaden
)
    map(Iterators.product(kpath, freq_list)) do (k, ω)
        momentum_set = momentum_calc(exciton, Q, k)
        k_h = momentum_set.k_h
        k_e = momentum_set.k_e
        k_i = momentum_set.k_i
        broaden(ω - E_exciton(exciton, Q) + E_v(exciton, k_h)) * abs2(wfn(k_i))
    end
end