using ProgressMeter
include("trion-solver.jl")
include("arpes-makie.jl")
using Printf

##

show_dispersion = true
σ = 10fs # Note that here σ tells us the width of the pulse; it should be *large* to produce δ-function like ARPES spectrum
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
kx_points = LinRange(-0.35, 1.7, 250)
ky_points = LinRange(-0.35, 1.7, 250)
dkx = step(kx_points)
dky = step(ky_points)
dk = dkx * dky
k_points = map(t -> SVector{2, Float64}(collect(t)), 
    collect(Iterators.product(kx_points, ky_points))
)
ω_points = LinRange(-1, 3, 500)
P_x = hexagonal_edge * 0.8
P_T  = SVector{2, Float64}([P_x, 0.0]) 
k_e = SVector{2, Float64}([-m_e / M * (w - P_x) + 0.1, 0.0])

ϵ_v2(k) =  - norm(k - w)^2 / 2m_h * inv_eV
ϵ_v1(k) =  - norm(k)^2     / 2m_h * inv_eV
ϵ_c1(k) =    norm(k)^2     / 2m_e * inv_eV + E_g
ϵ_c2(k) =    norm(k - w)^2 / 2m_e * inv_eV + E_g
E_SP(P) = inv_eV * norm(P .- w)^2 / 2M .+ E_B .+ 2E_g

M = 2m_h + m_e

ham = IndirectTwoBandModel2D(m_e, m_h, E_g, SVector{2, Float64}([w, 0.0]))
dielectric = Dielectric2D(ϵ)
broadening = gaussian_broadening(σ)

A_kω_Q = trionarpes_e1e2h1_thread(
    ham, E_B, 
    phi1sa1sb(IndirectTwoBandMat2D(ham, dielectric), a, b), 
    P_T, 
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


# The contributions of (various values of) kh1 and ω to the ARPES spectrum
# when ke is fixed to a given value.
Ake2ω_ke1 = map(Iterators.product(kx_points, ω_points)) do (k_e2x, ω)
    k_e1 = k_e
    k_e2 = SVector{2, Float64}([k_e2x, 0.0])
    k_h = P_T - k_e1 - k_e2
    ϵ_c2 =  inv_eV * (k_e2 - w)' * (k_e2 - w) / 2m_e
    ϵ_v  =  inv_eV * k_h' * k_h / 2m_h

    δ_factor = broadening(ω - E_S_PT - E_B - E_g + ϵ_c2 + ϵ_v)
            
    k_1 = k_e1 + m_e / M * w - m_e / M * P_T
    k_2 = k_e2 - (m_e + m_h) / M * w - m_e / M * P_T
    trion_structure_factor = abs(A_SQ_k1k2(k_1, k_2))^2
    
    δ_factor * trion_structure_factor
end

dispersion_Ake2ω_ke1 = map(kx_points) do k_e2x
    k_e1 = k_e
    k_e2 = SVector{2, Float64}([k_e2x, 0.0])
    k_h = P_T - k_e1 - k_e2
    ϵ_c2 =  inv_eV * (k_e2 - w)' * (k_e2 - w) / 2m_e
    ϵ_v  =  inv_eV * k_h' * k_h / 2m_h

    E_S_PT + E_B + E_g - ϵ_c2 - ϵ_v
end

Ake1ω_ke2 = map(Iterators.product(kx_points, ω_points)) do (k_e1x, ω)
    k_e2 = k_e
    k_e1 = SVector{2, Float64}([k_e1x, 0.0])
    k_h = P_T - k_e1 - k_e2
    ϵ_c1 =  inv_eV * k_e1' * k_e1 / 2m_e
    ϵ_v  =  inv_eV * k_h' * k_h / 2m_h

    δ_factor = broadening(ω - E_S_PT - E_B - E_g + ϵ_c1 + ϵ_v)
            
    k_1 = k_e1 + m_e / M * w - m_e / M * P_T
    k_2 = k_e2 - (m_e + m_h) / M * w - m_e / M * P_T
    trion_structure_factor = abs(A_SQ_k1k2(k_1, k_2))^2
    
    δ_factor * trion_structure_factor
end

dispersion_Ake1ω_ke2 = map(kx_points) do k_e1x
    k_e2 = k_e
    k_e1 = SVector{2, Float64}([k_e1x, 0.0])
    k_h = P_T - k_e1 - k_e2
    ϵ_c1 =  inv_eV * k_e1' * k_e1 / 2m_e
    ϵ_v  =  inv_eV * k_h' * k_h / 2m_h

    E_S_PT + E_B + E_g - ϵ_c1 - ϵ_v
end



##

fig = Figure(size=(1000,500))

ax_arpes_ω = Axis(fig[1, 1],
    ylabel = L"$ω$ (eV)",
    xlabel = L"A(k_\mathrm{e},\omega,k_\mathrm{e2})",
    xlabelsize = 20,
    ylabelsize = 20,
    xticklabelsize = 18,
    yticklabelsize = 18,
)

Aω_ke = sum(Ake2ω_ke1, dims=1)[1, :]
ω_Amax = ω_points[argmax(Aω_ke)]
hlines!(ax_arpes_ω, ω_Amax, color = colorant"gray80")
hlines!(ax_arpes_ω, maximum(dispersion_Ake2ω_ke1), color = colorant"gray80", linestyle = :dot)
lines!(ax_arpes_ω, Aω_ke, ω_points)

hidedecorations!(ax_arpes_ω, ticks = false, ticklabels = false, label = false)
ylims!(ax_arpes_ω, (minimum(ω_points), maximum(ω_points)))

ax_heatmap_single_ke2 = Axis(fig[1, 2], 
    xlabel = L"$k_{\mathrm{e2}}$ (Å)", 
    ylabel = L"$ω$ (eV)",
    xlabelsize = 20,
    ylabelsize = 20,
    xticklabelsize = 18,
    yticklabelsize = 18,
)

if show_dispersion
    lines!(ax_heatmap_single_ke2, kx_points, dispersion_Ake2ω_ke1, color = colorant"aqua")
end

hlines!(ax_heatmap_single_ke2, ω_Amax, color = colorant"gray80")
hlines!(ax_heatmap_single_ke2, maximum(dispersion_Ake2ω_ke1), color = colorant"gray80", linestyle = :dot)

heatmap!(ax_heatmap_single_ke2, kx_points, ω_points, Ake2ω_ke1,
    colormap = arpes_colormap(transparency_gradience),
)

ylims!(ax_heatmap_single_ke2, (minimum(ω_points), maximum(ω_points)))
hidedecorations!(ax_heatmap_single_ke2, ticks = false, ticklabels = false, label = false)

text!(ax_heatmap_single_ke2, 1.0, 0.2,
    fontsize = 14,
    text = latexstring(@sprintf "k_\\mathrm{e1} = %4.2f" k_e[1]), 
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
vlines!(ax_heatmap_single, -m_e / M * (w - P_T)[1], color = colorant"gray80")
vlines!(ax_heatmap_single, (m_e + m_h) / M * w[1] + m_e / M * P_x, color = colorant"gray80")

heatmap!(ax_heatmap_single, kx_points, ω_points, A_kω_Q,
    colormap = arpes_colormap(transparency_gradience),
)
hidedecorations!(ax_heatmap_single, ticks = false, ticklabels = false, label = false)
ylims!(ax_heatmap_single, (minimum(ω_points), maximum(ω_points)))

save("trion-arpes-eeh-ke1-$(@sprintf "%4.2f" k_e[1])-sigma-$(σ/fs)-nk-$nk_side.png", fig)

fig

##

let fig = Figure()
    Ake2_SQ_ke1 = map(k_points) do ke2
        ke1 = k_e
        k1 = ke1 + m_e / M * (w - P_T) 
        k2 = ke2 - (m_e + m_h) / M * w - m_e / M * P_T 
        A_SQ_k1k2(k1, k2)
    end
    ax = Axis(fig[1, 1], aspect=1)
    heatmap!(ax, kx_points, ky_points, Ake2_SQ_ke1)
    save("eeh-wfn.png", fig)
    fig
end

##

fig = Figure(size=(1000,500))

ax_arpes_ω = Axis(fig[1, 1],
    ylabel = L"$ω$ (eV)",
    xlabel = L"A(k_\mathrm{e},\omega,k_\mathrm{e1})",
    xlabelsize = 20,
    ylabelsize = 20,
    xticklabelsize = 18,
    yticklabelsize = 18,
)

Aω_ke = sum(Ake1ω_ke2, dims=1)[1, :]
ω_Amax = ω_points[argmax(Aω_ke)]
hlines!(ax_arpes_ω, ω_Amax, color = colorant"gray80")
hlines!(ax_arpes_ω, maximum(dispersion_Ake1ω_ke2), color = colorant"gray80", linestyle = :dot)
lines!(ax_arpes_ω, Aω_ke, ω_points)

hidedecorations!(ax_arpes_ω, ticks = false, ticklabels = false, label = false)
ylims!(ax_arpes_ω, (minimum(ω_points), maximum(ω_points)))

ax_heatmap_single_ke1 = Axis(fig[1, 2], 
    xlabel = L"$k_{\mathrm{e1}}$ (Å)", 
    ylabel = L"$ω$ (eV)",
    xlabelsize = 20,
    ylabelsize = 20,
    xticklabelsize = 18,
    yticklabelsize = 18,
)

if show_dispersion
    lines!(ax_heatmap_single_ke1, kx_points, dispersion_Ake1ω_ke2, color = colorant"aqua")
end

hlines!(ax_heatmap_single_ke1, ω_Amax, color = colorant"gray80")
hlines!(ax_heatmap_single_ke1, maximum(dispersion_Ake1ω_ke2), color = colorant"gray80", linestyle = :dot)

heatmap!(ax_heatmap_single_ke1, kx_points, ω_points, Ake1ω_ke2,
    colormap = arpes_colormap(transparency_gradience),
)

ylims!(ax_heatmap_single_ke1, (minimum(ω_points), maximum(ω_points)))
hidedecorations!(ax_heatmap_single_ke1, ticks = false, ticklabels = false, label = false)

text!(ax_heatmap_single_ke1, 1.0, 0.2,
    fontsize = 14,
    text = latexstring(@sprintf "k_\\mathrm{e2} = %4.2f" k_e[1]), 
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
vlines!(ax_heatmap_single, -m_e / M * (w - P_T)[1], color = colorant"gray80")
vlines!(ax_heatmap_single, (m_e + m_h) / M * w[1] + m_e / M * P_x, color = colorant"gray80")

heatmap!(ax_heatmap_single, kx_points, ω_points, A_kω_Q,
    colormap = arpes_colormap(transparency_gradience),
)
hidedecorations!(ax_heatmap_single, ticks = false, ticklabels = false, label = false)
ylims!(ax_heatmap_single, (minimum(ω_points), maximum(ω_points)))

save("trion-arpes-eeh-ke2-$(@sprintf "%4.2f" k_e[1])-Px-$(@sprintf "%4.2f" P_Tx)-sigma-$(σ/fs)-nk-$nk_side.png", fig)

fig

