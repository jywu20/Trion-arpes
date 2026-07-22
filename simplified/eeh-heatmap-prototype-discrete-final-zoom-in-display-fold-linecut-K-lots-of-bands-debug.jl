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
using Printf

caption_padding = 35

broaden = gaussian_broadening(50fs)
electron_color = :black
hole_color = :black

m_h = 0.21
m_e = 0.37
E_g = 2.67
w_side = 4π / (3 * 3.144817974)
w = SVector{2, Float64}([w_side, 0.0])
trion = Intervalley2DChandraTrion(
    m_h = m_h,
    m_e = m_e,
    w = w,
    E_g = E_g,
    # Binding energy for the ehh trion mode
    # should be E_B = 0.59 = 2.67 - 2.08 according to the band gap and the trion frequency from Nat. Comm. 2017 Dec 14;8(1):2117., 
    # however, here we're using Diana's exciton energy band data, 
    # the lowest exciton frequency at Q = 0 is 2.06 eV (run this script and the first line of console output is this number; note that this is the energy of the lowest optically active exciton mode, i.e. the A1s-like mode, which corresponds to the second lowest mode in the exciton band in PRL 115, 176801 (2015), because the lowest mode is a dark exciton mode),
    # and to maintain 2.13 - 2.08 = 0.05 eV (the former is again from Nat Comm. 2017 Dec 14;8(1):2117) difference between the lowest exciton mode and the trion mode in Nat Comm. 2017 Dec 14;8(1):2117,
    # we need to manually shift 0.59 upwards by 2.13 - 2.06 = 0.07 eV.
    E_B = 0.66,
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
shift = -0.0546272

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

println(minimum(eig_matrix_0[1, :]))
println(minimum(eig_matrix_K[1, :]))

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
    Avck_S ./= phase_factor
    Avck_S = Avck_S[1, 1, :] + Avck_S[1, 2, :]
    Avck_S ./ (Avck_S[1] / abs(Avck_S[1]))
end

#endregion 
##########################################

ω_list = LinRange(2.6, 2.6, 1)
kx_list = LinRange(w_side, w_side, 1)
k1_list = [SA[kx, 0.0] for kx in kx_list]

E_v1_curve = map(k1_list) do k_h
    -E_v1(trion, k_h)
end
E_c1_curve = map(k1_list) do k_e
    E_c1(trion, k_e)
end
E_v2_curve = map(k1_list) do k_h
    -E_v2(trion, k_h)
end
E_c2_curve = map(k1_list) do k_e
    E_c2(trion, k_e)
end

S_list_0 = [1, 2, 3, 4, 5, 6, 7, 8, ]
S_list_K = [1, 2, 3, 4,  ]
#S_list_0 = [1, 2, 3, 4, 5, 6,  ]
#S_list_K = [1, 2, 3, 4, 5, 6, ]
ik_K = argmin(map(i -> norm(rk[:, i]), 1 : size(rk)[2]))
S_list_0_spinor = []
S_list_K_spinor = []
for iS in 1 : 200 # The index is in Avck, i.e. my fully spinor calculation
    # Remove B excitons
    if norm(Avck[1, 1, :, iS] + Avck[1, 2, :, iS]) < 0.8
        continue
    end
    # Remove skin forbidden excitons
    if norm(Avck[1, 1, ik_K, iS]) < norm(Avck[1, 2, ik_K, iS])
        continue
    end

    append!(S_list_K_spinor, iS)
    # The exciton energies given in the eigenval_0_like files have a double degeneracy,
    # so we have to replicate the exciton wave function to match this.
    append!(S_list_0_spinor, iS)
    append!(S_list_0_spinor, iS)
    
    # First column: the index of the exciton mode after dark modes are removed.
    # Second column: the index of the exciton mode in the original Avck. 
    @printf "%3i  %3i  %8.6f   %8.6f   \n" length(S_list_0_spinor) iS eig_matrix_0[length(S_list_0_spinor), 7] (2E_g - trion.E_B - eig_matrix_0[length(S_list_0_spinor), 7])
end
S_list_0 = 1 : length(S_list_0_spinor)
S_list_K = 1 : length(S_list_K_spinor)


A1s_like = fetch_S(Avck, 2)
A2p_like_1 = fetch_S(Avck, 8)
A2p_like_2 = fetch_S(Avck, 6)
A2s_like = fetch_S(Avck, 10)

Ak1k2 = wfn(trion)
total = trion_ARPES_eeh(trion, P, Ak1k2, Homogeneous2DExciton, 
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
    (map(S_list_0_spinor) do iS
         # The sqrt(2) factor is due to the duplication of the Q=0 wave function mentioned above.
        fetch_S(Avck, iS) / sqrt(2)
    end)...,
    (map(S_list_K_spinor) do iS
        fetch_S(Avck, iS)
    end)...,
], [
    #rk, 
    #rk, 
    fill(rk, length(S_list_0))..., 
    fill(rk .+ [w_side, 0, 0], length(S_list_K))...,
], 
k1_list, ω_list, broaden)

let f = Figure()
    A_flat = fetch_S(Avck, 10)
    intensity_limit = maximum(abs.(A_flat))
    ax = Axis(f[1, 1], aspect=DataAspect())
    scatter!(ax, rk[1, :], rk[2, :], color=real.(A_flat), colorrange=(-intensity_limit, intensity_limit), colormap=cgrad(:balance),)
    f
end