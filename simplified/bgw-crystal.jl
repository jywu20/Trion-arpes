using HDF5
using LinearAlgebra

function read_b(name::String)
    fid = h5open(name)
    crystal = fid["mf_header/crystal/"]
    B = crystal["bvec"][] * crystal["blat"][]
    Bdot = crystal["bdot"][]
    @assert norm(B' * B - Bdot) < 1e-5
    B
end
