# PhotonicBandConnectivity.jl

Calculate minimum band connectivities below the fundamental gap in photonic crystals.

## Installation

You can install PhotonicBandConnectivity.jl from Julia's `pkg>` prompt (entered by typing `]` at the Julia REPL):

```jl
pkg> dev https://github.com/thchr/PhotonicBandConnectivity.jl
```

By `dev`ing it (rather than `add`ing), we also automatically install [SymmetryBases.jl](https://github.com/thchr/SymmetryBases.jl), which is a required (and, like this package, unregistered) dependency of PhotonicBandConnectivity.jl.

## Example usage

PhotonicBandConnectivity.jl provides `minimal_expansion_of_zero_freq_bands` that allows computing the minimum band connectivity below the fundamental of 3D photonic crystals. For example:
```jl
using PhotonicBandConnectivity

# choose a space group and whether time-reversal invariance exists
sgnum, has_tr = 86, true
# calculate minimum connectivity and symmetry decomposition
cⁱs, μᵀ, sb, idx¹ᴸ = minimal_expansion_of_zero_freq_bands(sgnum, timereversal=has_tr)
```

Which returns:
- `μᵀ`: the minimum band connectivity,
- `sb`: the Hilbert basis of the set of compatible band structures (via [SymmetryBases.jl](https://github.com/thchr/SymmetryBases.jl)),
- `cⁱs`: decomposition of the apolar solutions (nᴸ⁺ᵀ) in `sb`,
- `idx¹ᴸ`: index of the auxiliary longitudinal symmetry vector nᴸ in 'sb'.

These variables can be used to construct the transverse solutions' symmetry content via:

```jl
using Crystalline: prettyprint_symmetryvector

# unpack from indices in 'sb' to symmetry vectors:
nᴸ = sb[idx¹ᴸ]                          # longitudinal choice
nᴸ⁺ᵀs = [sum(sb[idxs]) for idxs in cⁱs] # associated apolar solutions

# build all transverse minimal connectivity solutions ('nᵀs'), then
# pick one ('nᵀ') and display it
nᵀs = [nᴸ⁺ᵀ - nᴸ for nᴸ⁺ᵀ in nᴸ⁺ᵀs]
nᵀ = nᵀs[1]
prettyprint_symmetryvector(stdout, nᵀ, sb.irlabs)
```

See additional examples in the [examples folder](https://github.com/thchr/PhotonicBandConnectivity.jl/tree/master/examples) and in the documentation strings of exported functions.


## Citation

If you find this package useful in your reseach, please cite our arXiv paper:

- T. Christensen, H.C. Po, J.D. Joannopoulos, & M. Soljačić, *Location and topology of the fundamental gap in photonic crystals*, [arXiv:2106.10267 (2021)](https://arxiv.org/abs/2106.10267).
