# The momentum unit is Å^-1. We are essentially dealing with wave vectors.
# Basically k∼2π/a.

# The mass unit is electron mass.

# The energy unit is eV. 
# This means the dispersion relation in our code should take the form of 
# E = ħ^2 k^2/2m * inv_eV (eV), where k is in Å^-1 and the unit of m is electron mass.
# Note that here the name "inv_eV" is from the fact that E is in eV,
# and the unit conversion constant should cancel with the eV unit: inv_eV * eV = 1.
# 
# So how to derive inv_eV? Basically:
# E in eV = ħ^2 k^2/(2m * eV) = ħ^2 (k in Å^-1)^2 Å^-2 / (2 (m in electron mass) * electron mass * eV)
# = ħ^2 Å^-2 / (electron mass * eV) * (k in Å^-1)^2 / 2 (m in electron mass)
# and therefore
# inv_eV = ħ^2 Å^-2 / (electron mass * eV) = 7.61996422.
const inv_eV = 7.61996422

# How many time units based on eV are there in one femtosecond.
# 1fs / (1eV / ħ) = 1.5193.
# Therefore, n fs = 1.5193n (eV / ħ).
const fs = 1.5193    