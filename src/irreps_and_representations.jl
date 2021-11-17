"""
    $(TYPEDSIGNATURES)

Get the compatibility basis `sb` associated with the `LGIrreps`s `lgirs` with or
without time-reversal symmetry (set by `timereversal=true` or `false`, respectively),
as well as the associated indexes into the Γ point irreps in `sb`, computed from the
Γ-point irreps `lgirs`.
The space group number is inferred from the provided vector of `LGIrrep`s.
"""
function compatibility_basis_and_Γidxs(lgirs::AbstractVector{LGIrrep{D}};
                                       timereversal::Bool=false, 
                                       allpaths::Bool=false) where D
    sgnum = num(group(first(lgirs)))
    # Find the Hilbert basis that respects the compatibility relations
    sb, _ = compatibility_basis(sgnum, D;
                        spinful=false, timereversal=timereversal, allpaths=allpaths)
    # Find the indices of the Γ irreps in `BRS::BandRepSet` and `sb::SymBasis` and how they
    # map to the corresponding irrep indices in `lgirs`.
    # TODO: note that the irrep-sorting in sb and lgirs is not always the same (e.g. in ±
    #       irreps), so we are not guaranteed that Γidxs is a simple range (e.g., it could 
    #       be [1,3,5,2,4,6]). We really ought to align the irreps sorting in `lgirreps`
    #       versus `bandreps` (BRS) and `compatibility_basis` (sb).
    Γidxs = get_Γidxs(lgirs, sb)

    return sb, Γidxs
end

function get_Γidxs(lgirs::AbstractVector{<:LGIrrep}, sb_or_brs::Union{BandRepSet, SymBasis})
    irlabs_sb_or_brs = irreplabels(sb_or_brs)
    irlabs_lgirs = Crystalline.formatirreplabel.(label.(lgirs))
    Γidxs = map(irlab->findfirst(==(irlab), irlabs_sb_or_brs), irlabs_lgirs)

    return Γidxs
end

# irrep-expansions/representation at Γ for the transverse (2T), longitudinal (1L), and triad
# (2T+1L) plane wave branches that touch ω=0 at Γ
"""
    find_representation²ᵀ⁺¹ᴸ(lgirs::AbstractVector{LGIrrep{3}})
    find_representation²ᵀ⁺¹ᴸ(sgnum::Integer; timereversal::Bool=true)
"""
function find_representation²ᵀ⁺¹ᴸ end
"""
    find_representation¹ᴸ(lgirs::AbstractVector{LGIrrep{3}})
    find_representation¹ᴸ(sgnum::Integer; timereversal::Bool=true)
"""
function find_representation¹ᴸ    end
"""
    find_representation²ᵀ(lgirs::AbstractVector{LGIrrep{3}})
    find_representation²ᵀ(sgnum::Integer; timereversal::Bool=true)
"""
function find_representation²ᵀ    end

for postfix in ("²ᵀ⁺¹ᴸ", "¹ᴸ", "²ᵀ")
    f = Symbol("find_representation"*postfix) # method to be defined
    symvals_fun = Symbol("get_symvals"*postfix)

    # "root" accessors via lgirs
    @eval function $f(lgirs::AbstractVector{<:Crystalline.AbstractIrrep{3}})
        lg = group(first(lgirs))
        symvals = $symvals_fun(lg)

        return find_representation(symvals, lgirs)
    end

    # convenience accessors via 
    @eval function $f(sgnum::Integer; timereversal::Bool=true)
        lgirs = lgirreps(sgnum, Val(3))["Γ"]
        timereversal && (lgirs = realify(lgirs))

        return $f(lgirs)
    end
end