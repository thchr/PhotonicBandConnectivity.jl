"""
    transverse_symmetry_vectors(sgnum, ::Val{D}[; timereversal::Bool = true, kws...])
    transverse_symmetry_vectors(sgnum, D[; timereversal::Bool = true, kws...])
    transverse_symmetry_vectors(brs, lgirsd[; separate_vrep::Bool = true])
                                                            --> Vector{SymmetryVector{D}}

Compute all "principal" solutions to the transverse photonic symmetry vector problem for
bands connected to ω=0. A principal solution is necessarily connected and cannot be
decomposed a sum of a Hilbert basis vector and another compatibility-respecting vector.

The principal solutions give all photonic bands that are intrinsically connected to ω=0.
The symmetry vector of any grouping of photonic bands can always be expressed as a sum
of an element from the principal solutions (or no element, if the group is not connected to
ω=0) and a non-negative integer span of the Hilbert basis vectors (see `compatibility_basis`
from SymmetryBases.jl).

The set of principal solutions include all minimum-connectivity solutions, but will also
include any higher-connectivity solutions which is intrinsically connected to ω=0. I.e.,
the principal solutions are a subset of the solutions obtained from 
[`minimal_expansion_of_zero_freq_bands`](@ref).

## Keyword arguments
- `timereversal` (default, `true`): specify whether there is time-reversal symmetry (`true`)
  or not (`false`); only applicable if the setting is specified via a space group number
  `sgnum` and dimensionality `D`.
- `separate_vrep` (default, `true`): if `true`, the singular virtual representation
  associated with the transverse modes at Γ and ω=0 will be introduced as a separate irrep,
  which is singly occupied (see [`transverse_vrep`](@ref)). The associated "ordinary" 
  constituent irrep multiplicites will be correspondingly subtracted from each
  `SymmetryVector`. In this case, all solutions will have strictly non-negative
  multiplicities.
"""
function transverse_symmetry_vectors(
    sgnum::Integer,
    Dᵛ::Val=Val(3);
    timereversal::Bool = true, 
    kws...)
    brs = calc_bandreps(sgnum, Dᵛ; timereversal)
    lgirsd = lgirreps(sgnum, Dᵛ)
    timereversal && realify!(lgirsd)
    return transverse_symmetry_vectors(brs, lgirsd; kws...)
end

function transverse_symmetry_vectors(sgnum::Integer, D::Integer; kws...)
    return transverse_symmetry_vectors(sgnum, Val(D); kws...)
end

function transverse_symmetry_vectors(
    brs::Union{BandRepSet, Collection{<:AbstractSymmetryVector{D}}},
    lgirsd::AbstractDict{<:AbstractString, Collection{LGIrrep{D}}};
    separate_vrep :: Bool = true
    ) where D

    F = smith(stack(brs))
    dᵇˢ = count(!iszero, F.SNF)
    S = @view F.S[:,1:dᵇˢ]

    # determine inhomogeneous constraints (nᵢ ≥ 0 for all kᵢ ≠ Γ and nᵢ ≥ n²ᵀΓᵢ for kᵢ = Γ)
    lgirsΓ = lgirsd["Γ"]
    n²ᵀΓ = find_representation²ᵀ(lgirsΓ) # NB irrep-sorting may differ from `brs.irlabs`
    inhom = zeros(Int, length(first(brs)))
    for (lgirᵢ, nᵢ) in zip(lgirsΓ, n²ᵀΓ)
        idx = something(findfirst(==(label(lgirᵢ)), irreplabels(brs)))
        inhom[idx] = nᵢ
    end

    # define translated cone
    tC_matrix = hcat(S, -inhom) # specifying `S*c ≥ inhom`
    # now solve for the lattice points in the inhomogeneous "cone", which is simply
    # a shifted cone (so not really a cone, but an unbounded polyhedron - but also a very
    # special kind of cone-like polyhedron; seems it's sometimes just called a "translated
    # cone". The homogeneous, untranslated cone is then the "recession" or "tail" cone of
    # the associated polyhedron
    C = PyNormaliz.Cone(inhom_inequalities = tC_matrix)
    C.Compute("HilbertBasis", "DualMode")
    lattice_points = C.LatticePoints()

    if !all(isone, @view lattice_points[:,end])
        error("unexpectedly obtained non-unity-graded lattice points: danger")
    end
    
    # go from the local "ℤ lattice" (into `S` columns) to the irrep-multiplicity lattice
    _ns = Ref(S) .* eachrow(@view lattice_points[:,1:end-1])
    ns = SymmetryVectors(_ns, irreplabels(brs), lgirsd)
    sort!(ns, by=occupation) # sort by occupation

    if separate_vrep
        # insert virtual rep explicitly, as a separate irrep, and then subtract the
        # associated "ordinary" constituent irrep multiplicites from each symmetry vector
        lgirsv = irreps(first(ns))
        Γidx = something(findfirst(lgirs->klabel(first(lgirs))=="Γ", lgirsv))
        lgirsΓ′ = lgirsv[Γidx]
        vrep = transverse_vrep(lgirsΓ′)
        n²ᵀΓ′ = find_representation²ᵀ(lgirsΓ′) # sorted as in `lgirsΓ′`, not `lgirsΓ`
        foreach(ns) do n
            multsΓ = n.multsv[Γidx]
            multsΓ .-= n²ᵀΓ′
            @assert all(≥(0), multsΓ)
            push!(n.multsv[Γidx], 1) # append virtual rep, singly occupied
        end
        push!(parent(lgirsv[Γidx]), vrep) # egal (===) for all `ns` ⇒ propagates to all `ns`
    end
    return ns
end

# ---------------------------------------------------------------------------------------- #

