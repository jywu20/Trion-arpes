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

broaden = gaussian_broadening(20fs)
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

files = ["finite-Q/eigenval_$(i)_likespin_plus_v_new" for i in 0:11]
Q_list = [
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

a = read_a("../../MoS2/MoS2/4-absorption-120-no-sym-gw/eigenvectors.h5") * au_in_angstrom
Q_length_list = map(Q_list) do Q
    Q[1] * 2π / a * sqrt(3)
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
all_eigs = [read_eigs(f) for f in files]

# Ensure all files have the same number of eigenvalues
nfiles = length(all_eigs)
neigs = length(all_eigs[1])
@assert all(length(e) == neigs for e in all_eigs) "Files have inconsistent number of eigenvalues!"

eig_matrix = hcat(all_eigs...) 

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
    Avck_S = Avck[:, :, :, iS]
    Avck_S = Avck_S[1, 1, :] + Avck_S[1, 2, :]
    reshape(Avck_S, (1, 1, length(Avck_S), 1))
end

#endregion 
##########################################

ω_list = LinRange(-0.3, 3.0, 200) 
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

S_list = [3, 4, 5, 6, 7, 8]

let f = Figure()
    Ak1k2 = wfn(trion)
    Akω_total = trion_ARPES_eeh(trion, P, Ak1k2, TwoBandTMDExciton, 
        [
            IntraValley2DExcitonHybridLow(exciton_direct),
            IntraValley2DExcitonHybridHigh(exciton_direct),
            (map(S_list) do iS
                Homogeneous2DExciton(Q_length_list, eig_matrix[4, eachindex(Q_length_list)])
            end)...,
            #Homogeneous2DExciton(Q_length_list, eig_matrix[4, eachindex(Q_length_list)]),
            exciton_K
        ], 
        # Note that we should NOT use the K momentum from the BGW run and convert it into Cartesian coordinates,
        # because it's in 1/au and not 1/Å. 
        [
            Avck_A1s_bright/2, 
            Avck_A1s_bright/2, 
            (map(S_list) do iS
                fetch_S(Avck, iS)
            end)...,
            Avck_A1s_bright,
        ], [
            rk, 
            rk, 
            fill(rk, length(S_list))..., 
            rk .+ [w_side, 0, 0],
        ], 
        k1_list, ω_list, broaden)

    ax = Axis(f[1, 1])
    colsize!(f.layout, 1, Aspect(1, 1))
    
    lines!(ax, kx_list, E_c1_curve, color=electron_color)
    lines!(ax, kx_list, E_c2_curve, color=electron_color)
    lines!(ax, kx_list, E_v1_curve, color=hole_color)
    lines!(ax, kx_list, E_v2_curve, color=hole_color)

    heatmap!(ax, kx_list, ω_list, Akω_total, colormap=arpes_colormap(transparency_gradience))
    ylims!(ax, (minimum(ω_list), maximum(ω_list)))
    hidedecorations!(ax, ticklabels = false, ticks = false)
    save("eeh-heatmap-prototype-discrete-final.png", f)

    f
end
