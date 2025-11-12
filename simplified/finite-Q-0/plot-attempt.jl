using CairoMakie
using LinearAlgebra
using HDF5
include("../exciton-band-gamma.jl")

# Find all files matching the naming pattern
files = ["eigenval_$(i)_likespin_plus_v_new" for i in 0:11]
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

function read_bdot(name::String)
    fid_ex = h5open(name)
    bmat = read(fid_ex["mf_header/crystal/blat"]) * read(fid_ex["mf_header/crystal/bvec"])
    bdot = read(fid_ex["mf_header/crystal/bdot"])
    close(fid_ex)
    @assert norm(bmat' * bmat - bdot) < 1e-10
    bdot[1:2, 1:2]
end

a = 3.144817974
Q_length_list = map(Q_list) do Q
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
all_eigs = [read_eigs(f) for f in files]

# Ensure all files have the same number of eigenvalues
nfiles = length(all_eigs)
neigs = length(all_eigs[1])
@assert all(length(e) == neigs for e in all_eigs) "Files have inconsistent number of eigenvalues!"

# Convert to a matrix: first index = file index, second = eigenvalue index
eig_matrix = hcat(all_eigs...) 

# Scissor shifts for convergence 
shift = -0.0546272
shift_2p_like = 0.0221814
shift_2s_like = 0.0828396

eig_matrix .+= shift
eig_matrix[3:6, :] .+= shift_2p_like
eig_matrix[7:end, :] .+= shift_2s_like

exciton = IntraValley2DExcitonHybrid(
    A = 0.6,
    m_e = 0.37,
    m_h = 0.21,
    α = -0.9,
    β = 4,
    E_g = 2.84,
    E_B = 2.84 - eig_matrix[1, 1],
)

θ = .5 * π

let first_n = 20
    nfiles = size(eig_matrix, 2)

    # Create a figure
    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1],
        xlabel = "Å⁻¹",
        ylabel = "Eigenvalue (eV)",
        title = "Eigenvalue evolution across files"
    )

    # Plot each eigenvalue line
    for j in 1:min(first_n, size(eig_matrix, 1))
        scatter!(ax, Q_length_list, eig_matrix[j, eachindex(Q_length_list)])
    end
    
    Qxs = LinRange(Q_length_list[1], Q_length_list[end], 100)
    lines!(ax, Qxs, map(Qxs) do Q
        E_exciton(IntraValley2DExcitonHybridLow(exciton), @SVector [Q, 0.0])
    end)
    lines!(ax, Qxs, map(Qxs) do Q
        E_exciton(IntraValley2DExcitonHybridHigh(exciton), @SVector [Q, 0.0])
    end)
    
    Qxs = LinRange(Q_length_list[1],0.2, 100)
    
    ex = Homogeneous2DExciton(Q_length_list, eig_matrix[1, eachindex(Q_length_list)])
    lines!(ax, Qxs, map(Q -> E_exciton(ex, @SVector [Q, 0.0]), Qxs))

    ex = Homogeneous2DExciton(Q_length_list, eig_matrix[2, eachindex(Q_length_list)])
    lines!(ax, Qxs, map(Q -> E_exciton(ex, @SVector [Q, 0.0]), Qxs))


    ex = Homogeneous2DExciton(Q_length_list, eig_matrix[3, eachindex(Q_length_list)])
    lines!(ax, Qxs, map(Q -> E_exciton(ex, @SVector [Q, 0.0]), Qxs))

    ex = Homogeneous2DExciton(Q_length_list, eig_matrix[10, eachindex(Q_length_list)])
    lines!(ax, Qxs, map(Q -> E_exciton(ex, @SVector [Q, 0.0]), Qxs))

    fig
end
