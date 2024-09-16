"""
$(TYPEDSIGNATURES)

Return whether the symmetry vector `nᵀ′` can represent an isolated grouping of transverse
photonic bands connected to zero-frequency. 

## Arguments

- `nᵀ′`: a symmetry vector of the ω=0 connected transverse bands. This vector should
  exclude the symmetry data of the singular ω=0 modes at Γ (e.g., if `nᵀ′` is a 2-band
  symmetry vector, its Γ-projection should be empty).
- `sb`: a compatibility Hilbert basis (see `compatibility_basis` in SymmetryBases).
- `lgirs`: a vector of `LGIrrep`s at Γ (see `lgirreps` in Crystalline).
- `F` or `brs`: the elementary band representations, provided either as a `BandRepSet` or as
  the Smith decomposition of its matrix form. Can be omitted (incurring then a slight extra
  cost).
"""
function is_transverse_bandstruct(
            nᵀ′::Vector{<:Integer},
            sb::SymBasis,
            lgirs::AbstractVector{LGIrrep{3}},
            F::Smith # Smith decomposition of EBR matrix
            )

    sb.compatbasis || error(DomainError(sb, "`sb` must be a basis for {BS}"))
    sb.spinful     && error(DomainError(sb, "`sb` must be a spinless basis"))

    nΓ¹ᴸ = find_representation¹ᴸ(lgirs)
    nΓ²ᵀ = find_representation²ᵀ(lgirs)
    
    Γidxs    = get_Γidxs(lgirs, sb)
    ntidxs¹ᴸ = find_symmetry_constrained_bases(sb, nΓ¹ᴸ, Γidxs)
    pick¹ᴸ   = argmin(fillings(sb)[ntidxs¹ᴸ])
    idx¹ᴸ    = ntidxs¹ᴸ[pick¹ᴸ]
    nᴸ       = sb[idx¹ᴸ]
    
    # we assume that `nᵀ′` refers to a transverse state _without_ the ω=0 irrep data at Γ,
    # so now we add back in the "surrogate" choice for Γ-irrep for the ω=0 modes
    nᵀ = copy(nᵀ′)
    nᵀ[Γidxs] .+= nΓ²ᵀ
    # then we build a possible regular band, by adding a longitudinal mode (that fulfills
    # the constraints on longitudinal modes at Γ)
    n = nᵀ + nᴸ

    # and finally, we test whether `n` is a band structure or not; if it is, then `n²ᵀ` is 
    # a valid transverse band structure
    return isbandstruct(n, F)
end

function is_transverse_bandstruct(
            nᵀ′::Vector{<:Integer},
            sb::SymBasis,
            lgirs::AbstractVector{LGIrrep{3}},
            brs::BandRepSet = bandreps(sb.sgnum, 3; 
                                       timereversal=sb.timeinvar, spinful=sb.spinful,
                                       allpaths=sb.allpaths)
            )

    return is_transverse_bandstruct(nᵀ′, sb, lgirs, smith(stack(brs)))
end