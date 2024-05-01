include("trion-solver.jl")
include("arpes-makie.jl")
using Printf

##

# Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
σ = 20fs
m_h = 0.21
m_e = 0.37
E_g = 2.84
ϵ = 6.4
E_B = -0.1
hexagonal_edge = 4π / (3 * 3.144817974)
w = hexagonal_edge
β = 1
a = 10.3
b = 25.2
nk_side = 150 
kx_points = LinRange(-0.7, 1.7, nk_side)
ky_points = LinRange(-0.7, 1.7, nk_side)
dkx = step(kx_points)
dky = step(ky_points)
dk = dkx * dky
k_points = map(t -> SVector{2, Float64}(collect(t)),
    collect(Iterators.product(kx_points, ky_points))
)
ω_points = LinRange(0, 3, 500)
P_Tx = 1.2 * w
P_T = SVector{2, Float64}([P_Tx, 0.0]) 

ϵ_v2(k) =  - norm(k - w)^2 / 2m_h * inv_eV
ϵ_v1(k) =  - norm(k)^2     / 2m_h * inv_eV
ϵ_c(k)  =    norm(k)^2     / 2m_e * inv_eV + E_g
E_SP(P) = inv_eV * norm(P .- w)^2 / 2M .+ E_B .+ E_g

M = 2m_h + m_e

ham = IndirectTwoBandModel2D(m_e, m_h, E_g, SVector{2, Float64}([w, 0.0]))
dielectric = Dielectric2D(ϵ)
broadening(x) = @fastmath exp(- σ^2 * x^2)

A_kω_Q = trionarpes_e1h1h2_thread(
    ham, E_B,
    phi1sa1sb(IndirectTwoBandMat2D(ham, dielectric), a, b), 
    SVector{2, Float64}([P_Tx, 0.0]), 
    k_points, ω_points, 
    map(kx -> SVector{2, Float64}([kx, 0.0]), 
        kx_points
    ),
    broadening
)     



m_e = ham.m_e
m_h = ham.m_h
E_g = ham.E_g
M = 2m_h + m_e
w = ham.w

E_S_PT = inv_eV * (P_T - w)' * (P_T - w) / 2M
A_SQ_k1k2 = phi1sa1sb(IndirectTwoBandMat2D(ham, dielectric), a, b)

k_e = SVector{2, Float64}([0.38, 0.0])

# The contributions of (various values of) kh1 and ω to the ARPES spectrum
# when ke is fixed to a given value.
Akh1ω_ke = map(Iterators.product(kx_points, ω_points)) do (k_h1x, ω) 
    k_h1 = SVector{2, Float64}([k_h1x, 0.0])
    k_h2 = P_T - k_e - k_h1
    ϵ_v_kh1 = - inv_eV * k_h1' * k_h1 / 2m_h
    ϵ_v_kh2 = - inv_eV * (k_h2 - w)' * (k_h2 - w) / 2m_h

    δ_factor = broadening(ω - E_S_PT - ϵ_v_kh1 - ϵ_v_kh2 - E_B - E_g)
    
    k_1 = k_h1 + m_h / M * w - m_h / M * P_T
    k_2 = k_h2 - (m_h + m_e) / M * w - m_h / M * P_T
    trion_structure_factor = abs(A_SQ_k1k2(k_1, k_2))^2
    
    δ_factor * trion_structure_factor
end

dispersion_Akh1ω_ke = map(kx_points) do k_h1x
    k_h1 = SVector{2, Float64}([k_h1x, 0.0])
    k_h2 = P_T - k_e - k_h1
    ϵ_v_kh1 = - inv_eV * k_h1' * k_h1 / 2m_h
    ϵ_v_kh2 = - inv_eV * (k_h2 - w)' * (k_h2 - w) / 2m_h

    E_S_PT + ϵ_v_kh1 + ϵ_v_kh2 + E_B + E_g
end

##

show_dispersion = true 

fig = Figure(size=(1000,500))

ax_arpes_ω = Axis(fig[1, 1],
    ylabel = L"$ω$ (eV)",
    xlabelsize = 20,
    ylabelsize = 20,
    xticklabelsize = 18,
    yticklabelsize = 18,
)

Aω_ke = sum(Akh1ω_ke, dims=1)[1, :]
ωke_Amax(ke) = - inv_eV * (P_Tx - w[1] - ke)^2 / 4m_h + E_B + E_g + E_S_PT
ω_Amax = ω_points[argmax(Aω_ke)]
hlines!(ax_arpes_ω, ω_Amax, color = colorant"gray80")
hlines!(
    ax_arpes_ω, 
    ωke_Amax(k_e[1]),
    color = colorant"gray80",
    linestyle = :dot 
)
lines!(ax_arpes_ω, Aω_ke, ω_points)

hidedecorations!(ax_arpes_ω, ticks = false, ticklabels = false, label = false)
ylims!(ax_arpes_ω, (minimum(ω_points), maximum(ω_points)))

ax_heatmap_single_kh1 = Axis(fig[1, 2], 
    xlabel = L"$k_{\mathrm{h1}}$ (Å)", 
    ylabel = L"$ω$ (eV)",
    xlabelsize = 20,
    ylabelsize = 20,
    xticklabelsize = 18,
    yticklabelsize = 18,
)


if show_dispersion
    lines!(ax_heatmap_single_kh1, kx_points, dispersion_Akh1ω_ke,
        color = colorant"aqua")
    vlines!(
        ax_heatmap_single_kh1, 
        (P_Tx - w[1] - k_e[1]) / 2, 
        color = colorant"gray80"
    )
    hlines!(
        ax_heatmap_single_kh1, 
        ωke_Amax(k_e[1]),
        color = colorant"gray80",
        linestyle = :dot 
    )
    hlines!(ax_heatmap_single_kh1, ω_Amax, color = colorant"gray80")
end

heatmap!(ax_heatmap_single_kh1, kx_points, ω_points, Akh1ω_ke,
    colormap = arpes_colormap(transparency_gradience),
)

let Ak1_SQ_ke = map(kx_points) do k1x
        k1 = SVector{2, Float64}([k1x, 0.0])
        k2 = m_e / M * (P_T - w) - k1 - k_e
        A_SQ_k1k2(k1, k2)
    end
    
    kh1x_ke = kx_points .+ m_h / M * (P_Tx - w[1])
    lines!(ax_heatmap_single_kh1, 
        kh1x_ke, 
        Ak1_SQ_ke * 20
    ) 

    sorted_indices = sortperm(Ak1_SQ_ke, rev=true)
    vlines!(ax_heatmap_single_kh1, kh1x_ke[sorted_indices[1:2]], color = colorant"gray80")
end

ylims!(ax_heatmap_single_kh1, (minimum(ω_points), maximum(ω_points)))
hidedecorations!(ax_heatmap_single_kh1, ticks = false, ticklabels = false, label = false)

text!(ax_heatmap_single_kh1, 1.0, 0.2,
    fontsize = 14,
    text = latexstring(@sprintf "k_\\mathrm{e} = %4.2f" k_e[1]), 
)

ax_heatmap_single = Axis(fig[1, 3], 
    xlabel = L"$k_{\mathrm{e}}$ (Å)", 
    ylabel = L"$ω$ (eV)",
    xlabelsize = 20,
    ylabelsize = 20,
    xticklabelsize = 18,
    yticklabelsize = 18,
)

vlines!(ax_heatmap_single, k_e[1], color = colorant"gray80")
hlines!(ax_heatmap_single, ω_Amax, color = colorant"gray80")
hlines!(
    ax_heatmap_single, 
    ωke_Amax(k_e[1]),
    color = colorant"gray80",
    linestyle = :dot 
)
lines!(ax_heatmap_single, kx_points, ωke_Amax.(kx_points), color = colorant"aqua")
heatmap!(ax_heatmap_single, kx_points, ω_points, A_kω_Q,
    colormap = arpes_colormap(transparency_gradience),
)
hidedecorations!(ax_heatmap_single, ticks = false, ticklabels = false, label = false)
ylims!(ax_heatmap_single, (minimum(ω_points), maximum(ω_points)))

save("trion-arpes-sigma-$(σ/fs)-nk-$nk_side.png", fig)

fig

##

let fig = Figure()
    idx_ke = argmin(abs.(kx_points .- k_e[1]))
    ax = Axis(fig[1, 1])
    lines!(ax, ω_points, A_kω_Q[idx_ke, :])
    lines!(ax, ω_points, Aω_ke)
    fig
end

##

let fig = Figure()
    Ak1_SQ_ke = map(k_points) do k1
        k2 = m_e / M * (P_T - w) - k1 - k_e
        A_SQ_k1k2(k1, k2)
    end
    ax = Axis(fig[1, 1])
    heatmap!(ax, kx_points, ky_points, Ak1_SQ_ke)
    fig
end