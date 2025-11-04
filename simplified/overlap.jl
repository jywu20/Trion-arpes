include("wfn.jl")
include("eeh.jl")
include("exciton.jl")
using HDF5

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

"""
Caveats:
- To avoid the problem of k path interpolation, the exciton modes given to this function 
  are to be represented by "valleys" in the exciton band structure:
  the structures of the exciton wave functions are assumed to be roughly the same in these valleys,
  and the only momentum resolved quantity is the exciton energy 
  (to be calculated by the `E_exciton(exciton, Q)` function).
- Following this, user should be aware of the valley structure of exciton wave functions supplied to the function:
  for instance, if we assume that the two electrons are in K and K' vallyes 
  and the hole is in the K valley,
  then we should NOT supply the K'-K' Q=0 exciton wave function to the calculation,
  or otherwise the function will mistaken the hole in K valley in the residue state of the trion
  and the hole in the K' valley in the exciton mode.
  That the holes in the exciton and in the trion residue state are in the same valley should be manually checked.
- It is suggested to use patched sampled exciton wave functions for better momentum spatial resolution.

The last dimension of `Avcks` is the exciton mode:
it's possible to have multiple exciton modes withe degenerate energies.
"""
function trion_ARPES_eeh(
    trion::Intervalley2DChandraTrion,
    P::SVector{2, Float64},
    wfn,
    ::Type{T},
    excitons::Vector{T},
    Avcks::Vector{Array{ComplexF64, 4}},
    rks::Vector{Matrix{Float64}},
    kpath::Vector{SVector{2, Float64}},
    freq_list::LinRange{Float64, Int64},
    broaden
) where {T <: TwoBandTMDExciton}
    tmap(Iterators.product(kpath, freq_list)) do (k, ω)
        sum(1 : size(Avcks)[1]) do iS
            rk = rks[iS]
            Avck = Avcks[iS]
            exciton = excitons[iS]
            # Exciton momentum. 
            # Note that due to BerkeleyGW's convention - see http://manual.berkeleygw.org/4.0/kernel-overview/ --
            # the total momentum of the exciton is -Q, and 
            # if the momentum of the conduction band electron is ke,
            # then the momentum of the valence band *electron* is ke + Q,
            # and thus the momentum of the valence band hole is - ke - Q.
            Q = P - k
            trans_mat = sum(1 : size(rk)[2]) do ik
                k_e_in_ex = SVector{2, Float64}(rk[1:2, ik])
                res = 0

                let momentum_set = momentum_calc_eeh_rel(trion, P, k, k_e_in_ex)
                    k_1 = momentum_set.k_1
                    k_2 = momentum_set.k_2

                    # There may be several exciton modes with degenerate energies
                    for iS′ in 1 : size(Avck)[4]
                        # [1, 1, ik, iS] means we're working on the lowest conduction band and highest valence band
                        res += Avck[1, 1, ik, iS′]' * wfn(k_1, k_2) 
                    end
                end
                
                let momentum_set = momentum_calc_eeh_rel(trion, P, k_e_in_ex, k)
                    k_1 = momentum_set.k_1
                    k_2 = momentum_set.k_2

                    for iS′ in 1 : size(Avck)[4]
                        res -= Avck[1, 1, ik, iS′]' * wfn(k_1, k_2) 
                    end
                end

                res
            end
            
            broaden(ω - E_trion_eeh(trion, P) + E_exciton(exciton, Q)) * abs2(trans_mat)
        end
    end
end