using HDF5
using Printf
using LinearAlgebra
using StaticArrays
include("units.jl")

function read_ex_wfc(name::String, k_shift::SVector{3, Float64})
    fid_ex = h5open(name)
    
    rk = fid_ex["mf_header/kpoints/rk"] |> read
    # [b1 b2 b3]
    bmat = read(fid_ex["mf_header/crystal/blat"]) * read(fid_ex["mf_header/crystal/bvec"])
    bdot = read(fid_ex["mf_header/crystal/bdot"])

    # Definition of metric
    @assert norm(bmat' * bmat - bdot) < 1e-10
    # Unit conversion: real space
    @assert abs(read(fid_ex["mf_header/crystal/alat"]) * au_in_angstrom - 3.144817974) < 1e-6
    # Unit conversion: relation between a and b
    @assert abs(read(fid_ex["mf_header/crystal/alat"]) - 2π / read(fid_ex["mf_header/crystal/blat"])) < 1e-7

    # Shift the origin of the coordinate system to the center of the momentum patch in question.
    ik_K = -1
    for i in 1 : size(rk)[2]
        if norm(rk[:, i] - k_shift) < 1e-4
            ik_K = i
            break
        end
    end
    k_K = rk[:, ik_K]
    rk .-= k_K
    # From crystal coordinates to Cartesian coordinates
    rk = bmat * rk
    # From atomic Rydberg units to Å
    rk /= au_in_angstrom
    
    Avck = fid_ex["exciton_data/eigenvectors"][1, 1, :, :, :, :, 1] + im*fid_ex["exciton_data/eigenvectors"][2, 1, :, :, :, :, 1] 
    # Normalization of wave function
    @assert abs(norm(Avck[:, :, :, 1]) - 1) < 1e-7
    
    (rk, Avck)
end



rk, Avck = read_ex_wfc("../../MoS2/MoS2/4-absorption-120-no-sym-gw/eigenvectors.h5", SVector{3, Float64}(0.333333333333333, 0.3333333333333, 0))
ik_K = argmin(map(i -> norm(rk[:, i]), 1 : size(rk)[2]))

for iS in 1 : 100
    # Remove B excitons
    if norm(Avck[1, 1, :, iS] + Avck[1, 2, :, iS]) < 0.8
        continue
    end
    # Remove skin forbidden excitons
    if norm(Avck[1, 1, ik_K, iS]) < norm(Avck[1, 2, ik_K, iS])
        continue
    end
    @printf "%3i   %6.4f   %6.4f  \n" iS norm(Avck[1, 1, :, iS]) norm(Avck[1, 2, :, iS])
end
