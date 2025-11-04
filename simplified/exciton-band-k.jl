include("units.jl")
include("model-data.jl")

@kwdef struct InterValley2DExcitonHybrid <: TwoBandTMDExciton
    m_e::Float64
    m_h::Float64

    α::Float64
    β::Float64
    
    E_g::Float64
    E_B::Float64

    w::SVector{2, Float64}
end

function E_exciton(ex::InterValley2DExcitonHybrid, Q::SVector{2, Float64})
    M = ex.m_e + ex.m_h
    (inv_eV / 2M + ex.α + ex.β) * norm(Q - ex.w)^2 + ex.E_g - ex.E_B
end
