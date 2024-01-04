include("../src/real-space.jl")
using Plots


let grid = UniformGrid2D(10, 10, 1000, 1000)
    ψ_1s = x -> ϕ_r_1s_def(x, 1.5)
    overlap_numerical(ψ_1s, ψ_1s, grid)
end

let Ns = [50, 100, 200, 500, 1000, 2000, 4000, 5000],
    Ls = [10, 20, 30],
    ψ_1s = x -> ϕ_r_1s_def(x, 2.8),
    p = plot()
    for L in Ls
        overlap = map(Ns) do N
            grid = UniformGrid2D(L, L, N, N)
            overlap_numerical(ψ_1s, ψ_1s, grid)
        end
        plot!(p, log.(Ns), overlap, label = L, markershape = :circle)
    end
    p # It's expected that the values go to 1 when the parameters are large
end

let Ns = [50, 100, 200, 500, 1000, 2000, 4000, 5000],
    Ls = [10, 20, 30],
    a = 2.0,
    ψ_1s = x -> ϕ_r_1s_def(x, a),
    p = plot()
    for L in Ls
        E_coulomb = map(Ns) do N
            grid = UniformGrid2D(L, L, N, N)
            coulomb_2D_numerical(ψ_1s, ψ_1s, grid)
        end
        plot!(p, log.(Ns), E_coulomb, label = L, markershape = :circle)
    end
    hline!(p, [2 / a], label = "Analytic")
    p # It's expected that the values go to 2/a when the parameters are large
end

let Ls = [10, 20, 30],
    Ns = [50, 100, 200],    
    a = 0.9,
    ψ_1s = x -> ϕ_r_1s_def(x, a)
    p = plot()
    
    for L in Ls
        E_coulomb = map(Ns) do N
            grid = UniformGrid2D(L, L, N, N)
            coulomb_twopoint_2D_numerical(ψ_1s, ψ_1s, grid)
        end
        plot!(p, log.(Ns), E_coulomb, label = L, markershape = :circle)
    end
    p
end