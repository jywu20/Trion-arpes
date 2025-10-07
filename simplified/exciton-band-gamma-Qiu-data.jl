include("exciton-band-gamma.jl")
using CairoMakie

exciton = IntraValley2DExcitonHybrid(
    A = 0.6,
    m_e = 0.37,
    m_h = 0.21,
    α = -0.9,
    β = 4,
    E_g = 2.84,
    E_B = 0.71,
)

θ = .5 * π

Qxs = LinRange(-0.1, 0.1, 100)
Qs = map(Qxs) do Qx
    @SVector [Qx, 0.0]
end
E1 = map(Qs) do Q
    solve(exciton, Q, θ).values[1]
end
E2 = map(Qs) do Q
    solve(exciton, Q, θ).values[2]
end
P1 = map(Qs) do Q
    abs2(solve(exciton, Q, θ).vectors[1, 1])
end
P2 = map(Qs) do Q
    abs2(solve(exciton, Q, θ).vectors[1, 2])
end


let f = Figure()
    ax = Axis(f[1, 1])
    scatter!(ax, Qxs, E1)
    scatter!(ax, Qxs, E2)
    f
end

##

let f = Figure()
    ax = Axis(f[1, 1])
    scatter!(ax, Qxs, P1)
    ylims!(ax, (0, 1))
    f
end


##

let f = Figure()
    ax = Axis(f[1, 1])
    scatter!(ax, Qxs, map(Qs) do Q
        E_exciton(IntraValley2DExcitonHybridLow(exciton), Q)
    end)
    scatter!(ax, Qxs, map(Qs) do Q
        E_exciton(IntraValley2DExcitonHybridHigh(exciton), Q)
    end)
    save("exciton-bands.png", f)
    f
end