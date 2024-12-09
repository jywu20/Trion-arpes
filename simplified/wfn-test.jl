include("wfn.jl")
using Test
using ProgressMeter

@testset "Numerical stability of @fastmath in 1s wave function" begin
    a_ex_list = [0.5, 1.0, 2.0]
    kx_list = LinRange(0, 5, 10)
    progress = Progress(length(a_ex_list) * length(kx_list)^2)
    
    for a_ex in a_ex_list
        for kx in kx_list, ky in kx_list
            k = SVector{2, Float64}([kx, ky])
            @test ϕ_1s(k, a_ex) ≈ ϕ_1s_def(k, a_ex)
            next!(progress)
        end
    end
end
