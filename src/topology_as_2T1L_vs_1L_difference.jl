"""
    calc_topology_singular(
            nᵀ⁺ᴸ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer}, 
            BRS_B_F::Union{BandRepSet, AbstractMatrix{<:Integer}, Smith})
                                                                     -> ::TopologyKind

Determines whether a transverse symmetry vector (``T``) - defined as the difference of a
transverse + longitudinal (``T+L``) symmetry vector `nᵀ⁺ᴸ` and a longitudinal (``L``)
symmetry vector `nᴸ` - is topologically trivial or nontrivial from a symmetry perspective.

**Input:**
- `nᵀ⁺ᴸ`: can be computed from an index vector `cⁱ` and a compatibility Hilbert basis
  `sb::SymBasis` by `sum_symbases(sb, cⁱ)`;
- `nᴸ = sb[idx¹ᴸ]`: where `idx¹ᴸ` is the 1L-constraints pick in the `sb` basis;
- `BRS_B_F`: a matrix representation of the elementary band representations (EBR) basis,
  [typically obtained from `bandreps(sgnum, kwarg...)`]. Must be of type `BandRepSet`,
  `AbstractMatrix{<:Integer}`, or `Smith` (in order of decreasing conversion-related
  overhead).

**Output:**
- a member of the enum `SymBases.TopologyKind`, either `TRIVIAL = 0` or `NONTRIVIAL = 1`.

# Implementation
A transverse, singular band is trivial if [`indicators_singular`](@ref) returns a trivial set
of topological indices ``ν``, i.e., if ``νᵢ = 0`` for all ``i``. Otherwise, it is stably
nontrivial.
The indices computed from `indicators_singular` are well-defined and robust, because while
neither the ``T+L`` or ``L`` symmetry vectors in general are unique, their difference is
(i.e. ``T``, up to singular ``Γ``-irrep content): and as is the difference of their
topological indices.
"""
function calc_topology_singular(
            nᵀ⁺ᴸ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer}, F::Smith)
    inds, _ = indicators_singular(nᵀ⁺ᴸ, nᴸ, F)

    return all(iszero, inds) ? TRIVIAL : NONTRIVIAL
end
@deprecate topology_from_2T1L_xor_1L calc_topology_singular

# ---------------------------------------------------------------------------------------- #

"""
    indicators_singular(
            nᵀ⁺ᴸ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer}, F::Smith)
    --> Vector{Int}, AbstractVector{Int}

Returns the indices and the nontrivial factors corresponding to `nᵀ⁺ᴸ`, `nᴸ`, and `F` for a
zero-frequency (singular) solution.
"""
function indicators_singular(
            nᵀ⁺ᴸ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer}, F::Smith)

    # get topological indices of L and T+L symmetry vectors
    indicesᴸ, Λ   = indicators(nᴸ,   F)
    indicesᵀ⁺ᴸ, _ = indicators(nᵀ⁺ᴸ, F)

    # define 2T solution's topology from Λ-modded difference:
    indicesᵀ = mod.(indicesᵀ⁺ᴸ .- indicesᴸ, Λ)

    return indicesᵀ, Λ
end

# ---------------------------------------------------------------------------------------- #
# CONVENIENCE ACCESSORS/WRAPPERS

for f in (:calc_topology_singular, :indicators_singular)
    # Basic wrappers
    @eval function $f(nᵀ⁺ᴸ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer},
                      B::AbstractMatrix{<:Integer})
        F = smith(B)            # ::Smith
        return $f(nᵀ⁺ᴸ, nᴸ, F)
    end

    @eval function $f(nᵀ⁺ᴸ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer},
                      BRS::BandRepSet)
        B = matrix(BRS)         # ::Matrix{Int}
        return $f(nᵀ⁺ᴸ, nᴸ, B)
    end

    # Fancier, more convenient wrappers
    @eval begin 
        @doc """
            $($f)(nᵀ::Vector{<:Integer}, nᴸ::Vector{<:Integer}, 
                  m²ᵀ::AbstractVector{<:Integer}, Γidxs::AbstractVector{<:Integer}, 
                  BRS_B_F::Union{BandRepSet, AbstractMatrix{<:Integer}, Smith})
                                                                           -> ::TopologyKind
        
        A convenience wrapper over
        [`$($f)(::AbstractVector{<:Integer}, ::AbstractVector{<:Integer}, ::Smith)`](@ref),
        taking the _physical_ (i.e. _without_ singular Γ-irrep content) transverse symmetry
        vector `nᵀ`, the pinned 2T Γ-symmetry content `nΓ²ᵀ`, and the Γ-irrep indexing
        vector `Γidx` instead of `nᵀ⁺ᴸ` from which it reconstructs the latter (using the
        provided longitudinal symmetry `nᴸ` as well). This is frequently a more convenient
        access point.
        
        As for its 3-argument parent method, the final argument `BRS_B_F` must provide the
        EBR basis as a `BandRepSet`, an `AbstractMatrix{<:Integer}`, or a `Smith`
        decomposition.
        """
        function $f(nᵀ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer},
                    nΓ²ᵀ::AbstractVector{<:Integer}, Γidxs::AbstractVector{<:Integer},
                    BRS_B_F::Union{BandRepSet, AbstractMatrix{<:Integer}, Smith})

            nᵀ⁺ᴸ = copy(nᵀ)      # reconstruct nᵀ⁺ᴸ from nᵀ, nᴸ, and nΓ²ᵀ (note the tricky
            nᵀ⁺ᴸ[Γidxs] .+= nΓ²ᵀ # broadcasted indexing with Γidxs)
            nᵀ⁺ᴸ .+= nᴸ

            return $f(nᵀ⁺ᴸ, nᴸ, BRS_B_F)
        end
    end

    # TODO: doc-string (convenience accessor that doesn't require us to actually provide 
    #       `nᴸ` or `ms²ᵀ`)
    @eval function $f(nᵀ::AbstractVector{<:Integer},
                      sb::SymBasis,
                      lgirs::AbstractVector{LGIrrep{3}}, # Γ-irreps
                      BRS_B_F::Union{BandRepSet, AbstractMatrix{<:Integer}, Smith})

        sb.compatbasis || error(DomainError(sb, "`sb` must be a basis for {BS}"))
        sb.spinful     && error(DomainError(sb, "`sb` must be a spinless basis"))

        nΓ¹ᴸ = find_representation¹ᴸ(lgirs)
        nΓ²ᵀ = find_representation²ᵀ(lgirs)

        Γidxs    = get_Γidxs(lgirs, sb)
        ntidxs¹ᴸ = find_symmetry_constrained_bases(sb, nΓ¹ᴸ, Γidxs)
        pick¹ᴸ   = argmin(fillings(sb)[ntidxs¹ᴸ])
        idx¹ᴸ    = ntidxs¹ᴸ[pick¹ᴸ]
        nᴸ       = sb[idx¹ᴸ]

        return $f(nᵀ, nᴸ, nΓ²ᵀ, Γidxs, BRS_B_F)
    end
end
