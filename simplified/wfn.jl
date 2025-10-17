# Wave function of the internal degrees of freedom of the trion.

using StaticArrays
using LinearAlgebra
include("model-data.jl")

# The line below previously was
# ϕ_1s_def(k::SVector{2, Float64}, a_ex::Float64) = 1 / (1 + norm(k)^2 * a_ex^2 / 4)^2
# This is possibly a previously unreported bug.
# There are two problems. 
#
# First, the exponent should be 3/2 for a 2D system:
# The 3D wave function takes the form of ∫ r^2 sin θ dr dθ dφ e^(-r/a),
# which gives us an exponent of 2.
# On the other hand, for a 2D system, what we have is ∫ r dr dθ e^(-r/a), 
# and using Mathematica we can easily verify that the exponent is actually 3/2.
#
# The second problem is, it should be 1 + k^2 a^2, and not 1 + k^2 a^2 / 4.
# Our ansatz should be consistent with that in  Theory of neutral and charged excitons in monolayer transition metal dichalcogenides,
# where Eq. (5) reads 2/πa^2 exp(−ρ/a), and not 2ρ/a.
# It takes a form like exp(-2ρ/a) for Fourier transform to produce 
# (4 + k^2 a^2)^1.5, and hence 1 + k^2 a^2 / 4.
ϕ_1s_def(k::SVector{2, Float64}, a_ex::Float64) = 1 / (1 + norm(k)^2 * a_ex^2)^1.5
@inline ϕ_1s(k::SVector{2, Float64}, a_ex::Float64) = 
    @fastmath 1 / (1 + k' * k * a_ex^2)^1.5

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

function wfn(exciton::InterValley2DExciton)
    a = exciton.a
    @inline ϕ_1s_a(k) = ϕ_1s(k, a)
    ϕ_1s_a
end