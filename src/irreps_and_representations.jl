"""
    $(TYPEDSIGNATURES)

Get the compatibility basis `sb` associated with the `LGIrreps`s `lgirs` with or
without time-reversal symmetry (set by `timereversal=true` or `false`, respectively),
as well as the associated indexes into the Γ point irreps in `sb`, computed from the
Γ-point irreps `lgirs`.
The space group number is inferred from the provided vector of `LGIrrep`s.
"""
function compatibility_basis_and_Γidxs(
    lgirs::AbstractVector{LGIrrep{D}};
    timereversal::Bool=false, 
    allpaths::Bool=false
) where D
    sgnum = num(group(first(lgirs)))
    # Find the Hilbert basis that respects the compatibility relations
    sb, _ = compatibility_basis(sgnum, D;
                        spinful=false, timereversal=timereversal, allpaths=allpaths)
    # Find the indices of the Γ irreps in `sb::SymBasis` and how they map to the
    # corresponding irrep indices in `lgirs`.
    # TODO: note that the irrep-sorting in sb and lgirs is not always the same (e.g. in ±
    #       irreps), so we are not guaranteed that Γidxs is a simple range (e.g., it could 
    #       be [1,3,5,2,4,6]). We really ought to align the irreps sorting in `lgirreps`
    #       versus `bandreps` (brs) and `compatibility_basis` (sb).
    Γidxs = get_Γidxs(lgirs, sb)

    return sb, Γidxs
end

function get_Γidxs(
    lgirs::AbstractVector{LGIrrep{D}},
    sb_or_brs::Union{Collection{SpinlessBandRep{D}}, SymBasis{D}}
) where D
    irlabs_sb_or_brs = irreplabels(sb_or_brs)
    irlabs_lgirs = label.(lgirs)
    Γidxs = map(irlab->findfirst(==(irlab), irlabs_sb_or_brs), irlabs_lgirs)

    return Γidxs
end

# irrep-expansions/representation at Γ for the transverse (2T), longitudinal (1L), and triad
# (2T+1L) plane wave branches that touch ω=0 at Γ
"""
    find_representation²ᵀ⁺¹ᴸ(lgirs::AbstractVector{LGIrrep{D}}, symval_optargs...)
    find_representation²ᵀ⁺¹ᴸ(sgnum::Integer; timereversal::Bool=true, symval_optargs...)
"""
function find_representation²ᵀ⁺¹ᴸ end
"""
    find_representation¹ᴸ(lgirs::AbstractVector{LGIrrep{D}})
    find_representation¹ᴸ(sgnum::Integer; timereversal::Bool=true)
"""
function find_representation¹ᴸ    end
"""
    find_representation²ᵀ(lgirs::AbstractVector{LGIrrep{D}}, symval_optargs...)
    find_representation²ᵀ(sgnum::Integer; timereversal::Bool=true, symval_optargs...)
"""
function find_representation²ᵀ    end

for postfix in ("²ᵀ⁺¹ᴸ", "¹ᴸ", "²ᵀ")
    f = Symbol("find_representation"*postfix) # method to be defined
    symvals_fun = Symbol("get_symvals"*postfix)

    # "root" accessors via lgirs
    @eval function $f(
        lgirs::AbstractVector{<:Crystalline.AbstractIrrep{D}},
        symval_optargs...
    ) where D
        lg = group(first(lgirs))
        symvals = $symvals_fun(lg, symval_optargs...)

        return find_representation(symvals, lgirs)
    end

    # convenience accessors via 
    @eval function $f(
        sgnum::Integer,
        ::Val{D}=Val(3);
        timereversal::Bool=true,
        symval_optargs...
    ) where D
        lgirs = lgirreps(sgnum, Val(D))["Γ"]
        timereversal && (lgirs = realify(lgirs))

        return $f(lgirs, symval_optargs...)
    end
end