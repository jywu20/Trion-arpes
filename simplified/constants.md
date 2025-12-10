Constants used in the calculation
===========

# Quasiparticle band structure

Parameters about the quasiparticle band structure are from [PRL 111,216805 (2013)](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.111.216805).
From Table I, we find:
- the effective mass of the electron: 0.37 electron mass
- the effective mass of the hole: 0.21 electron mass
- band gap:  2.84 eV

Note that the 2.84 eV band gap is not fully converged;
it should be 2.67 eV according to [a newer paper](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.93.235435).

# Comparison of (unconverged) trion and exciton BSE in Nat. Comm.  2117 (2017)

Parameters about trion binding are from [Nat. Comm.  2117 (2017)](https://www.nature.com/articles/s41467-017-02286-6).
In Table 1, we find:
- The fundamental band gap is 2.89 eV (slightly larger than the band gap above but it's ok)
- The negative trion frequency is 2.08 eV
- The positive trion frequency is 2.09 eV
- The exciton frequency is 2.13 eV

These parameters listed above are calculated without a substrate.
But the parameters with substrate are not too different.

# Full convergence of exciton BSE in PRL 111,216805 (2013)

There is probably a convergence issue:
the lowest exciton energy reported in [PRL 111,216805 (2013)](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.111.216805)
is 1.88 eV (see Table II).
The largest k-grid used in [Nat. Comm.  2117 (2017)](https://www.nature.com/articles/s41467-017-02286-6) is 39 39 1 (see Supplementary Fig. 3);
the largest k-grid used in [PRL 111,216805 (2013)](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.111.216805) is 72 72 1,
so the latter is more reliable.

# Comparison between direct and K-K' excitons in PRL 115, 176801

In [PRL 115, 176801](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.176801),
it is mentioned that the direct band gap is 2.67 eV
and the binding energy of the unlike-spin transition state,
which is the lowest, is 0.65 eV (in the main text).
Furthermore, the lowest energy Q = K exciton is 5 meV lower in energy than the lowest Q = 0 exciton (see the supplementary material, below Fig. 1).

# Variational ansatz of excitons and trions

The trion variational parameters are from [PRB 88,045318 (2013)](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.88.045318).
They are:  a = 10.3 Å and b = 25.2 Å.

The exciton variational parameters is from [PRB 88,045318 (2013)](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.88.045318). It is 10.4 Å.

Note that the binding energy is defined differently in [PRB 88,045318 (2013)](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.88.045318),
where we find that the binding energy is reported to be 26 meV (see Table 1).
This data comes from [Nat. Mat. 12, pages207–211 (2013)](https://www.nature.com/articles/nmat3505),
where the trion binding energy $E_{\text{A}^-}$ is defined in Eq. (1):
$$
\omega_{\text{A}} - \omega_{\text{A}^-} = E_{\text{A}^-} + E_{\text{F}}.
$$
You can see that the binding energy is just the difference
between the exciton A mode and the trion A- mode,
and we need to remove the Fermi energy. 