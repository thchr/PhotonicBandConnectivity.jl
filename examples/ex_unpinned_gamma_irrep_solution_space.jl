# Example of how to determine the physically admissible solution space of the zero-frequency
# Γ-irreps of 3D photonic crystals (i.e., computing rows 4, 5, and 6 of Tables S7-S8 in 
# [Phys. Rev. X **12**, 021066](https://doi.org/10.1103/PhysRevX.12.021066))

# ---------------------------------------------------------------------------------------- #
using PhotonicBandConnectivity
using Crystalline

sgnum = 215
lgirs = realify(lgirreps(sgnum)["Γ"])
nfixed, nfree = physical_zero_frequency_gamma_irreps(
    lgirs; 
    supergroup_constraints=true,
    force_fixed=true,
    lattice_reduce=true)
prettyprint_irrep_solspace(stdout, nfixed, nfree, lgirs)