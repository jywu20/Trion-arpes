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

ϕ_1s_def(k::SVector{2, Float64}, a_ex::Float64) = 1 / (1 + norm(k)^2 * a_ex^2 / 4)^2
@inline ϕ_1s(k::SVector{2, Float64}, a_ex::Float64) = @fastmath 1 / (1 + k' * k * a_ex^2 / 4)^2