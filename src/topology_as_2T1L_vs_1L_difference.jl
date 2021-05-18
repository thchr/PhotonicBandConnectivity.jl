"""
    topology_from_2T1L_xor_1L(
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
- a member of the enum `SymBases.TopologyKind`, either `trivial=0` or `nontrivial=1`.

# Implementation
A simple "xor" approach is followed: the topology of the``T`` solutions are inferred from
the "difference" of the topologies of the ``T+L`` and ``L`` solutions in the xor (`⊻`)
sense, i.e. ``T`` is nontrivial if the topologies of ``T+L`` and ``L`` differ:

    ┌────────────┬────────────╥────────────╖        ┌───┬───╥───────╖
    │     L      │    T+L     ║    T Z₂    ║        │ x │ y ║ x ⊻ y ║
    ├────────────┼────────────╫────────────╢        ├───┼───╫───────╢
    |  trivial   │  trivial   ║  trivial   ║        | 0 │ 0 ║   0   ║
    |  trivial   │ nontrivial ║ nontrivial ║        | 0 │ 1 ║   1   ║
    | nontrivial │  trivial   ║ nontrivial ║        | 1 │ 0 ║   1   ║
    | nontrivial │ nontrivial ║  trivial   ║        | 1 │ 1 ║   0   ║
    └────────────┴────────────╨────────────╜        └───┴───╨───────╜

This is a meaningful definition, because while neither the ``T+L`` or ``L`` symmetry vectors
in general are unique, their difference is (i.e. ``T``, up to singular ``Γ``-irrep content):
same goes for the associated topology.
"""
function topology_from_2T1L_xor_1L(
            nᵀ⁺ᴸ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer}, F::Smith)
    # check topology of L and T+L symmetry vectors
    is_nontrivialᴸ   = calc_topology(nᴸ,   F) == NONTRIVIAL
    is_nontrivialᵀ⁺ᴸ = calc_topology(nᵀ⁺ᴸ, F) == NONTRIVIAL
    # infer 2T solution's topology from xor-difference:
    is_nontrivialᵀ = is_nontrivialᴸ ⊻ is_nontrivialᵀ⁺ᴸ

    return is_nontrivialᵀ ? NONTRIVIAL : TRIVIAL
end

# convenience accessors
function topology_from_2T1L_xor_1L(
            nᵀ⁺ᴸ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer},
            B::AbstractMatrix{<:Integer})
    F = smith(B)
    return topology_from_2T1L_xor_1L(nᵀ⁺ᴸ, nᴸ, F)
end

function topology_from_2T1L_xor_1L(
            nᵀ⁺ᴸ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer},
            BRS::BandRepSet)
    B = matrix(BRS, true)                   # ::Matrix{Int}
    return topology_from_2T1L_xor_1L(nᵀ⁺ᴸ, nᴸ, B)
end

# -----------------------------------------------------------------------------------------

"""
    topology_from_2T1L_xor_1L(nᵀ::Vector{<:Integer}, nᴸ::Vector{<:Integer}, 
                m²ᵀ::AbstractVector{<:Integer}, Γidxs::AbstractVector{<:Integer}, 
                BRS_B_F::Union{BandRepSet, AbstractMatrix{<:Integer}, Smith})
                                                                -> ::TopologyKind

A wrapper over the equivalently named 3-argument signature method, but taking the _physical_
(i.e. _without_ singular Γ-irrep content) transverse symmetry vector `nᵀ`, the pinned 2T
symmetry content `m²ᵀ`, and the Γ indexing vector `Γidx` instead of `nᵀ⁺ᴸ` from which it
reconstructs the latter.
This can often be a more convenient access point.

In similarity with its 3-argument parent method, its final argument `BRS_B_F` must be a
representation of the elementary band representation (EBR) basis, given as either a
`BandRepSet`, an `AbstractMatrix{<:Integer}`, or a `Smith` decomposition.
"""
function topology_from_2T1L_xor_1L(
            nᵀ::AbstractVector{<:Integer}, nᴸ::AbstractVector{<:Integer},
            m²ᵀ::AbstractVector{<:Integer}, Γidxs::AbstractVector{<:Integer},
            BRS_B_F::Union{BandRepSet, AbstractMatrix{<:Integer}, Smith})

    nᵀ⁺ᴸ = copy(nᵀ)      # reconstruct T+L solution from nᵀ, nᴸ, and m²ᵀ (note the tricky
    nᵀ⁺ᴸ[Γidxs] .+= m²ᵀ  # broadcasted indexing with Γidxs)
    nᵀ⁺ᴸ .+= nᴸ

    return topology_from_2T1L_xor_1L(nᵀ⁺ᴸ, nᴸ, BRS_B_F)
end

# TODO: doc-string (convenience accessor that doesn't require us to actually provide `nᴸ` or
#       `ms²ᵀ`)
function topology_from_2T1L_xor_1L(
            nᵀ::AbstractVector{<:Integer},
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
    nᴸ       = sb[pick¹ᴸ]

    topology_from_2T1L_xor_1L(nᵀ, nᴸ, nΓ²ᵀ, Γidxs, BRS_B_F)
end