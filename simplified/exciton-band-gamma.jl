using StaticArrays
using Dierckx
using LinearAlgebra
include("units.jl")
include("model-data.jl")

@kwdef struct IntraValley2DExcitonHybrid
    A::Float64

    m_e::Float64
    m_h::Float64

    α::Float64
    β::Float64
    
    E_g::Float64
    E_B::Float64
end

struct IntraValley2DExcitonHybridLow <: TwoBandTMDExciton
    data::IntraValley2DExcitonHybrid
end

struct IntraValley2DExcitonHybridHigh <: TwoBandTMDExciton
    data::IntraValley2DExcitonHybrid
end

σ0 = @SMatrix [1 0; 0 1]
σx = @SMatrix [0 1; 1 0]
σy = @SMatrix [0 -im; im 0]

function solve(exciton::IntraValley2DExcitonHybrid, Q::SVector{2, Float64}, θ::Float64)
    A = exciton.A
    M = exciton.m_e + exciton.m_h
    α = exciton.α
    β = exciton.β
    σ = cos(2θ) * σx + sin(2θ) * σy
    H = A * (σ0 +  σ) * norm(Q)
    H += ((inv_eV / 2M + α + β) * σ0 + β * σ) * norm(Q, 2)^2
    eigen(H)
end

function E_exciton(ex::IntraValley2DExcitonHybridLow, Q::SVector{2, Float64})
    solve(ex.data, Q, 0.0).values[1] + ex.data.E_g - ex.data.E_B
end

function E_exciton(ex::IntraValley2DExcitonHybridHigh, Q::SVector{2, Float64})
    solve(ex.data, Q, 0.0).values[2] + ex.data.E_g - ex.data.E_B
end

struct Homogeneous2DExciton <: TwoBandTMDExciton
    Q_norm::Vector{Float64}
    band::Vector{Float64}
    extrapolate::Spline1D
end

function Homogeneous2DExciton(Q_norm, band)
    Homogeneous2DExciton(Q_norm, band, Spline1D(Q_norm, band, bc="extrapolate"))
end

function E_exciton(ex::Homogeneous2DExciton, Q::SVector{2, Float64})
    Q_norm = norm(Q)
    ex.extrapolate(Q_norm)
end
