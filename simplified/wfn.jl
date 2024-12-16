# Wave function of the internal degrees of freedom of the trion.

using StaticArrays
using LinearAlgebra
include("model-data.jl")

ϕ_1s_def(k::SVector{2, Float64}, a_ex::Float64) = 1 / (1 + norm(k)^2 * a_ex^2 / 4)^2

@inline ϕ_1s(k::SVector{2, Float64}, a_ex::Float64) = 
    @fastmath 1 / (1 + k' * k * a_ex^2 / 4)^2

function wfn(trion::Intervalley2DChandraTrion) 
    a = trion.a
    b = trion.b

    @inline ϕ_1s_a(k) = ϕ_1s(k, a)  
    @inline ϕ_1s_b(k) = ϕ_1s(k, b)
    @inline A_SQ_k1k2(k_1, k_2) = 1/sqrt(2) * (
        ϕ_1s_a(k_1) * ϕ_1s_b(k_2) + ϕ_1s_a(k_2) * ϕ_1s_b(k_1)
    )
    
    A_SQ_k1k2
end

function wfn(exciton::IntraValley2DExciton)
    a = exciton.a
    @inline ϕ_1s_a(k) = ϕ_1s(k, a)
    ϕ_1s_a
end