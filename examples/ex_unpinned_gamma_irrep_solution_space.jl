# Example of how to determine the physically admissible solution space of the zero-frequency
# Γ-irreps of 3D photonic crystals (i.e., computing rows 4, 5, and 6 of Tables S7-S8 in 
# [Phys. Rev. X **12**, 021066](https://doi.org/10.1103/PhysRevX.12.021066))
# See also `physical_zero_frequency_gamma_irreps_O3` which is faster (but whose
# implementation relies on the insights of the thinking behind the implementation of
# `physical_zero_frequency_gamma_irreps`).
# ---------------------------------------------------------------------------------------- #
using PhotonicBandConnectivity
using Crystalline

sgnum = 215
lgirs = realify(lgirreps(sgnum)["Γ"])
nfixed, nfree = PhotonicBandConnectivity.physical_zero_frequency_gamma_irreps(
    lgirs; 
    supergroup_constraints=true,
    force_fixed=true,
    lattice_reduce=true)
prettyprint_irrep_solspace(stdout, nfixed, nfree, lgirs)