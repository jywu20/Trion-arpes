# Plot the ARPES heatmap;

include("units.jl")
include("qp-bands.jl")
include("broadening.jl")
include("exciton.jl")
include("ehh.jl")
include("wfn.jl")
include("arpes.jl")
include("overlap.jl")
include("exciton-band-gamma.jl")
include("exciton-band-k.jl")
include("bgw-crystal.jl")
using CairoMakie
using LaTeXStrings
using Colors

caption_padding = 35

broaden = gaussian_broadening(40fs)
electron_color = colorant"deepskyblue2"
hole_color = colorant"coral2"

m_h = 0.21
m_e = 0.37
E_g = 2.84
w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])
trion = Intervalley2DChandraTrion(
    m_h = m_h,
    m_e = m_e,
    w = w,
    E_g = E_g,
    E_B = 0.75, # Binding energy for the ehh trion mode
    a = 10.3,
    b = 25.2
)

##########################################
#region Reading BerkeleyGW outputs 

k_K = SVector{3, Float64}(0.333333333333333, 0.3333333333333, 0)
rk, Avck = read_ex_wfc("../../MoS2/MoS2/4-absorption-120-no-sym-gw/eigenvectors.h5", k_K)
k_K_real = read_B("../../MoS2/MoS2/4-absorption-120-no-sym-gw/eigenvectors.h5") * k_K

files_0 = ["finite-Q-0/eigenval_$(i)_likespin_plus_v_new" for i in 0:11]
files_K = ["finite-Q-K/eigenval_$(i)_unlikespin_plus_v_new" for i in 0:13]
Q_list_0 = [
    [0.0,0.0],
    [0.0005,0.0005],
    [0.001,0.001],
    [0.003,0.003],
    [0.005,0.005],
    [0.0075,0.0075],
    [0.01,0.01],
    [0.012,0.012],
    [0.015,0.015],
    [0.02,0.02]
]

Q_list_K = [
    [-0.02, -0.02],
    [-0.015, -0.015],
    [-0.01, -0.01],
    [-0.005, -0.005],
    [-0.003, -0.003],
    [-0.001, -0.001],
    [0.001, 0.001],
    [0.003, 0.003],
    [0.005, 0.005],
    [0.0075, 0.0075],
    [0.01, 0.01],
    [0.012, 0.012],
    [0.015, 0.015],
    [0.02, 0.02]
]

a = read_a("../../MoS2/MoS2/4-absorption-120-no-sym-gw/eigenvectors.h5") * au_in_angstrom
Q_length_list_0 = map(Q_list_0) do Q
    (2 / sqrt(3)) * Q[1] * 2π / a * sqrt(3)
end
Q_length_list_K = map(Q_list_K) do Q
    (2 / sqrt(3)) * Q[1] * 2π / a * sqrt(3) 
end


# Read eigenvalues from a file
function read_eigs(filename)
    eigs = Float64[]
    open(filename, "r") do io
        for line in eachline(io)
            # Skip comment lines
            startswith(line, "#") && continue
            fields = split(strip(line))
            isempty(fields) && continue
            # Parse the first number (eigenvalue column)
            push!(eigs, parse(Float64, fields[1]))
        end
    end
    return eigs
end

# Read all eigenvalue lists
all_eigs_0 = [read_eigs(f) for f in files_0]
all_eigs_K = [read_eigs(f) for f in files_K]

# Ensure all files have the same number of eigenvalues
nfiles = length(all_eigs_0)
neigs = length(all_eigs_0[1])
@assert all(length(e) == neigs for e in all_eigs_0) "Files have inconsistent number of eigenvalues!"

eig_matrix_0 = hcat(all_eigs_0...) 
eig_matrix_K = hcat(all_eigs_K...) 

# Scissor shifts for convergence 
#shift = -0.0546272
shift = 0 # The trion binding energy we use here should be lower than the exciton binding energy, 
# so we have to use slightly unconverged exciton energies.
# This can be solved by obtaining a better estimation of the trion energy
# but I won't bother.
shift_2p_like = 0.0221814
shift_2s_like = 0.0828396

eig_matrix_0 .+= shift
eig_matrix_0[3:6, :] .+= shift_2p_like
eig_matrix_0[7:end, :] .+= shift_2s_like

shift_2p_unlike = 0.0233762
shift_2s_unlike = 0.0746934

eig_matrix_K .+= shift
eig_matrix_K[2:3, :] .+= shift_2p_unlike
eig_matrix_K[2:3, :] .-= (eig_matrix_K[2:3, :] .- eig_matrix_K[2, 1]) * 0.2
eig_matrix_K[4:end, :] .+= shift_2s_unlike

#endregion
############################################

##########################################
#region The linear-parabolic bands

exciton_direct = IntraValley2DExcitonHybrid(
    A = 0.6,
    m_e = 0.37,
    m_h = 0.21,
    α = -0.9,
    β = 4,
    E_g = E_g,
    E_B = 0.71,
)

exciton_K = InterValley2DExcitonHybrid(
    m_e = 0.37,
    m_h = 0.21,
    α = -0.9,
    β = 4,
    E_g = E_g,
    E_B = 0.71,
    w = w
)

# The trion momentum is set to be w
P_ratio = 1.
P = P_ratio * w


# The second lowest mode is bright for K; we need to concentrate on that mode.
Avck_A1s_bright = Avck[:, :, :, 2]
# In Fig. 1bcd in PRL 115, 176801 (2015), 
# we note that the spin up and spin down bands cross with each other:
# in a patch around K, the lowest conduction band is spin up,
# but around it, the lowest conduction band is spin down,
# while the highest valence band is spin up. 
# In our trion two-band model, we expect the conduction and valence bands to have the same spin,
# and this means the conduction band in our two-band model corresponds to 
# c1 near the patch around K and c2 around that patch.
# 
# Therefore, after reading the exciton wave function, 
# we have to reorder the bands. 
# But there is one way to cheat:
# because the lowest optically active exciton mode is quite accurately two-band,
# we can just do the follows:
Avck_A1s_bright = Avck_A1s_bright[1, 1, :] + Avck_A1s_bright[1, 2, :]
Avck_A1s_bright = reshape(Avck_A1s_bright, (1, 1, length(Avck_A1s_bright), 1))
println(size(Avck_A1s_bright))

function fetch_S(Avck::Array{ComplexF64, 4}, iS::Int)
    phase_factor = Avck[:, :, :, 2] ./ abs.(Avck[:, :, :, 2])
    Avck_S = Avck[:, :, :, iS]
    Avck_S = Avck_S[1, 1, :] + Avck_S[1, 2, :]
    A_vck_S = reshape(Avck_S, (1, 1, length(Avck_S), 1)) ./ phase_factor
    A_vck_S ./ (A_vck_S[1] / abs(Avck_S[1]))
end

#endregion 
##########################################

ω_list = LinRange(2.0, 3.0, 200) 
kx_list = LinRange(-0.35, 1.7, 250) 
k1_list = [SA[kx, 0.0] for kx in kx_list]

E_c1_curve = map(k1_list) do k_h
    -E_v1(trion, k_h)
end
E_v1_curve = map(k1_list) do k_e
    E_c1(trion, k_e)
end
E_c2_curve = map(k1_list) do k_h
    -E_v2(trion, k_h)
end
E_v2_curve = map(k1_list) do k_e
    E_c2(trion, k_e)
end

S_list_0 = [1, 2, 3, 4, 5, 6, 7, 8, ]
S_list_K = [1, 2, 3, 4,  ]
#S_list_0 = [1, 2, 3, 4, 5, 6,  ]
#S_list_K = [1, 2, 3, 4, 5, 6, ]

A1s_like = fetch_S(Avck, 2)
A2p_like_1 = fetch_S(Avck, 8)
A2p_like_2 = fetch_S(Avck, 6)
A2s_like = fetch_S(Avck, 10)

let f = Figure(size=(600, 400))
    Ak1k2 = wfn(trion)
    Akω_total = trion_ARPES_eeh(trion, P, Ak1k2, Homogeneous2DExciton, 
        [
            #IntraValley2DExcitonHybridLow(exciton_direct),
            #IntraValley2DExcitonHybridHigh(exciton_direct),
            (map(S_list_0) do iS
                Homogeneous2DExciton(Q_length_list_0, eig_matrix_0[iS, eachindex(Q_length_list_0)])
            end)...,
            #Homogeneous2DExciton(Q_length_list_0, eig_matrix_0[4, eachindex(Q_length_list_0)]),
            #exciton_K
            map(S_list_K) do iS
                Homogeneous2DExciton(Q_length_list_K, eig_matrix_K[iS, eachindex(Q_length_list_K)], shift=w)
            end...,
        ], 
        # Note that we should NOT use the K momentum from the BGW run and convert it into Cartesian coordinates,
        # because it's in 1/au and not 1/Å. 
        [
            # There should be a 1/2 factor for the first two wave functions,
            # because the lowest K and K' excitons are hybridized and the form of the resulting wave function 
            # has been analytically found in https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.115.176801.
            # Still, we expect similar hybridization to happen in higher states,
            # where we don't really know the coefficients.
            # Therefore we just omit the 1/2 factor to avoid introducing non physical intensity difference
            # between the low and high excitons.
            (map([2, 2, 6, 6, 8, 8, 10, 10]) do iS
                fetch_S(Avck, iS) / sqrt(2)
            end)...,
            (map([2, 6, 8, 10]) do iS
                fetch_S(Avck, iS)
            end)...,
        ], [
            #rk, 
            #rk, 
            fill(rk, length(S_list_0))..., 
            fill(rk .+ [w_side, 0, 0], length(S_list_K))...,
        ], 
        k1_list, ω_list, broaden)

    ax = Axis(f[1, 1])
    #colsize!(f.layout, 1, Aspect(1, 1))
    
    lines!(ax, kx_list, E_c1_curve, color=electron_color)
    lines!(ax, kx_list, E_c2_curve, color=electron_color)
    lines!(ax, kx_list, E_v1_curve, color=hole_color)
    lines!(ax, kx_list, E_v2_curve, color=hole_color)

    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))

    kx_list_0 = LinRange(w_side-0.08, w_side+0.08, 100)
    for S in S_list_0
        ex = Homogeneous2DExciton(Q_length_list_0, eig_matrix_0[S, eachindex(Q_length_list_0)])
        lines!(ax, kx_list_0, map(kx_list_0) do kx
            E_trion_eeh(trion, P) - E_exciton(ex, (P - @SVector [kx, 0]))
        end, color = colorant"darkorange", linestyle=:dash)
    end
    kx_list_K = LinRange(-0.08, 0.08, 100)
    for S in S_list_K
        ex = Homogeneous2DExciton(Q_length_list_K, eig_matrix_K[S, eachindex(Q_length_list_K)], shift=w)
        lines!(ax, kx_list_K, map(kx_list_K) do kx
            E_trion_eeh(trion, P) - E_exciton(ex, (P - @SVector [kx, 0]))
        end, color = colorant"darkorange", linestyle=:dash)
    end

    ylims!(ax, (minimum(ω_list), maximum(ω_list)))
    hidedecorations!(ax, ticklabels = false, ticks = false)
    save("eeh-heatmap-prototype-discrete-final-zoom-in.png", f)

    f
end
