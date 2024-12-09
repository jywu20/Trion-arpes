using StaticArrays

"""
Two-band description of trion properties of 2D valleytronics  materials in which 
- There are two valley-peak pairs, where the spin degrees of freedom are strongly locked to the valley degrees of freedom (and there is no non-trivial spin texture)
- The positive and negative trion modes are all intervalley, which means the spin part (which is also the valley part) of the trion wave function is a singlet state.

The effective masses at the valley and the peak are `m_e` and `m_h`, respectively.
Note that here we assume that `m_h` is positive,
and the negative signs should be manually inserted when necessary.
`w` gives the distance between the two peaks in the momentum space.

When the dielectric properties of the system is given, the trion properties can be calculated.
A good approximation is Chandrasekhar's ansatz.
Currently we don't have the relevant solvers,
so the two radii are to be manually picked.
"""
@kwdef struct Intervalley2DChandraTrion
    # The two-band model
    m_e::Float64
    m_h::Float64
    w::SVector{2, Float64}
    E_g::Float64

    # Info about the trion modes;
    # with Chandrasekhar's ansatz, the wave function of the positive and negative trions
    # are the same 
    E_B::Float64
    a::Float64
    b::Float64
end

