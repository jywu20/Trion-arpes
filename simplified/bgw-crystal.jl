using HDF5
using LinearAlgebra

function read_B(name::String)
    fid = h5open(name)
    crystal = fid["mf_header/crystal/"]
    B = crystal["bvec"][] * crystal["blat"][]
    Bdot = crystal["bdot"][]
    @assert norm(B' * B - Bdot) < 1e-5
    B
end

function read_a(name::String)
    fid_ex = h5open(name)
    bmat = read(fid_ex["mf_header/crystal/blat"]) * read(fid_ex["mf_header/crystal/bvec"])
    bdot = read(fid_ex["mf_header/crystal/bdot"])
    @assert norm(bmat' * bmat - bdot) < 1e-10
    a = read(fid_ex["mf_header/crystal/alat"])
    close(fid_ex)
    a
end