using Crystalline
using Crystalline: AbstractIrrep, find_parent_pointgroup
using LinearAlgebra
using LLLplus

# Implements `physical_zero_frequency_gamma_irreps` which returns a solution space for the
# physically admissible zero-frequency Γ-irreps, given a set of Γ-irreps.

# ---------------------------------------------------------------------------------------- #

"""
    physical_zero_frequency_gamma_irreps(
        lgirs                  :: AbstractVector{LGIrrep{3}};
        supergroup_constraints :: Bool = true,
        force_fixed            :: Bool = true,
        lattice_reduce         :: Bool = true
        )                                       --> Tuple{Vector{Int}, Matrix{Int}}

Compute and return a solution space for the physically admissible zero-frequency Γ-irreps,
given a set of Γ-irreps `lgirs`.

The solution space is returned as a tuple `(nfixed :: Vector{Int}, nfree :: Matrix{Int})`,
parametrizing the admissible solutions according to:
    
        `nfixed + nfree * [l, p, q, …]`

where `[l, p, q, …]` is a vector of free integer variables. The rows of `nfixed` and `nfree`
correpond to multiplicities of the irreps in `lgirs`.

The solution space can be pretty-printed with [`prettyprint_irrep_solspace`](@ref).

# Implementation limitations
The inputs `lgirs` must be a set of `realified` or inherently real irreps; time-reversal
breaking (i.e., complex or pseudoreal) irreps are not yet supported.

# Keyword arguments
- `supergroup_constraints :: Bool = true`: if `true`, imposes the constraints not just from
  the parent point group of `lgirs` little group, but also from its maximal super point 
  group.  constrained by
  `supergroup_constraints = true` (`false`) is equivalent to the "strong" ("weak") 
  constraints listed in Tables S7-S8 of 
  [Phys. Rev. X **12**, 021066](https://doi.org/10.1103/PhysRevX.12.021066).
- `force_fixed :: Bool = true`: if `true`, the fixed part of the solution space (`nfixed`)
  is chosen to match that from `PhotonicBandConnectivity.find_representation²ᵀ(lgirs))`; 
  i.e., the fixed part agrees with column 4 of Tables S7-S8.
- `lattice_reduce :: Bool = true`: if `true`, attempt to simplify the choice of basis terms
  in the free term `nfree` (i.e., the lattice term of the solution space) using
  lattice-reduction techniques (Seysen algorithm).
"""
function physical_zero_frequency_gamma_irreps(
    lgirs                  :: AbstractVector{LGIrrep{3}};
    supergroup_constraints :: Bool = true,
    force_fixed            :: Bool = true,
    lattice_reduce         :: Bool = true,
)
    timereversal = true
    if any(lgir->reality(lgir) ≠ REAL && !lgir.iscorep, lgirs)
        # COMPLEX or PSEUDOREAL irreps in `lgirs`, but neither converted to corep form;
        # i.e., we must have `timereversal = false`
        timereversal = false
    end
    lg = group(first(lgirs))
    if klabel(lg) != "Γ" # only Γ-irreps are supported
        error(lazy"input `lgirs` are associated with the $(label(lg))-point: must be Γ")
    end
    pg, Iᵖ²ᵍ, _ = find_isomorphic_parent_pointgroup(lg)
    pgirs = pgirreps(label(pg))
    timereversal && (pgirs = realify(pgirs))
    # TODO: `timereversal = false` is actually still not yet supported, since it will fail
    #       later, when we assume the character table consists of real integers.
    pgops = operations(group(first(pgirs)))
    if Iᵖ²ᵍ ≠ 1:length(lg)
        # ops are permuted between `pg` & `lg`; re-sort `pg` & `pgirs` (we intentionally do
        # not use `permute!` below, because the data in `pgirs` is unsafely loaded via JLD2
        # and mutation causes problems)
        pgops = pgops[Iᵖ²ᵍ]
        pg′ = PointGroup{3}(pg.num, pg.label, pgops)
        for (i, pgir) in enumerate(pgirs)
            matrices′ = pgir.matrices[Iᵖ²ᵍ]
            pgirs[i] = PGIrrep{3}(pgir.cdml, pg′, matrices′, pgir.reality, pgir.iscorep)
        end
    end

    local fixed, free, unpinned_idxmap
    if supergroup_constraints
        spg = find_pg_supergroup(pg)
        spgirs = pgirreps(label(spg), Val(3))
        timereversal && (spgirs = realify(spgirs))
        spgops = operations(group(first(spgirs)))
        # NB: ops-permutation between `spg`, `spgops` & `lg` is no problem; we re-sort below
        
        # find unpinned operations & indexes into them in `ops`
        _, sis = pinned_and_unpinned_operation_indexes(spgops)
        # find fixed and free parts of possible symeig solutions
        fixed, free, _ = unpinned_symeigs_basis(spgirs, lattice_reduce)
        # compute an index map between operations in the referenced in the rows of `free` &
        # `fixed`, versus the operations in `pgops` (and hence, `pgirs`)
        unpinned_idxmap = [findfirst(==(op), @views spgops[sis]) for op in pgops]
    else
        _, is = pinned_and_unpinned_operation_indexes(pgops)
        fixed, free, _ = unpinned_symeigs_basis(pgirs, lattice_reduce)
        unpinned_idxmap = [findfirst(==(op), @views pgops[is]) for op in pgops]
    end

    # prepare to find ambiguity in terms of Γ-irreps rather than symmetry eigenvalues
    x = Vector{Int}(undef, length(pgops))
    xfree = zeros(Int, length(pgops), size(free, 2))
    for (i, op) in enumerate(pgops)
        if !isnothing(unpinned_idxmap[i]) # operation with unpinned symeig
            i′ = unpinned_idxmap[i]
            x[i] = fixed[i′]
            xfree[i,:] = free[i′,:]
        else                              # operation with pinned symeig
            x[i] = pinned_symval²ᵀ(op)
        end
    end

    # force the fixed part to be that chosen by PhotonicBandConnectivity
    if force_fixed
        x_force = round.(Int, PhotonicBandConnectivity.get_symvals²ᵀ(pgops))
        # verify that there are indeed integer solutions also with `x_force`
        xfreeᵍ, maxλ_xfreeᵍ = generalized_inv(xfree)
        @assert all(isinteger, (xfreeᵍ*(x_force - x))./maxλ_xfreeᵍ)
        x = x_force
    end

    # now we have the symeigs lattice that we want; convert it to an irrep multipl. lattice
    ct = characters(pgirs)
    nfixed = ct\x
    nfree = ct\xfree
    @assert isapprox(nfixed, round.(real(nfixed)), atol=1e-10)  
    nfixed = round.(Int, real(nfixed))
    @assert isapprox(nfree, round.(real(nfree)), atol=1e-10) 
    nfree = round.(Int, real(nfree))

    # check that irrep labels have the same sorting input irreps (`lgirs`) and point group
    # irreps (`pgirs`); the rows of `nfixed` and `nfree` are presently in the sorting of
    # `pgirs`, so may potentially need a permutation if the sorting is different
    pgirlabs = ct.irlabs
    lgirlabs = label.(lgirs)
    if pgirlabs == lgirlabs
        return nfixed, nfree
    else
        # we need to find a permutation of the irrep labels between `lgirs` and `pgirs` &
        # permute `nΓ` and `nΓfree` accordingly
        if Set(pgirlabs) != Set(lgirlabs)
            error("irrep labels of `lgirs` and `pgirs` differ by more than a permutation")
        end
        perm = [something(findfirst(==(lgirlab), pgirlabs)) for lgirlab in lgirlabs]
        return nfixed[perm], nfree[perm, :]
    end
end

function physical_zero_frequency_gamma_irreps(sgnum :: Integer; kws...)
    physical_zero_frequency_gamma_irreps(realify(lgirreps(sgnum, Val(3))["Γ"]); kws...)
end

# ---------------------------------------------------------------------------------------- #
# Pretty-printing

"""
    prettyprint_irrep_solspace(
        [io::IO,]
        nfixed::AbstractVector,
        nfree::AbstractMatrix,
        lgirs::Union{AbstractVector{<:LGIrrep}, AbstractVector{String}}
        )

Pretty-print the solution space obtained by `physical_zero_frequency_gamma_irreps` to `io`.
The solution space is parametrized as:

    `nfixed + nfree * [l, p, q, …]`

where `[l, p, q, …]` is a vector of free integer variables.
"""
function prettyprint_irrep_solspace(
            io::IO,
            nfixed::AbstractVector,
            nfree::AbstractMatrix,
            irlabs::AbstractVector{String})
    size(nfree, 2) > 3 && error("Unexpectedly high number of free variables; unhandled")
    vars = ('l','p','q')
    Crystalline.prettyprint_symmetryvector(io, nfixed, irlabs; braces=false)
    for (i,n) in enumerate(eachcol(nfree))
        iszero(n) && continue
        print(io, " + ")
        print(io, "(")
        Crystalline.prettyprint_symmetryvector(io, n, irlabs; braces=false)
        print(io, ")", vars[i])
    end
end
function prettyprint_irrep_solspace(
            io::IO,
            nfixed::AbstractVector,
            nfree::AbstractMatrix,
            lgirs::AbstractVector{<:LGIrrep})
    prettyprint_irrep_solspace(io, nfixed, nfree, label.(lgirs))
end

# ---------------------------------------------------------------------------------------- #

function unpinned_symeigs_basis(
            irs::AbstractVector{<:AbstractIrrep{3}},
            lattice_reduce :: Bool = true
            )
    ct = characters(irs)
    ops = operations(ct)
    χ = ct
    χ′ = round.(Int, real(ct)) # NB/TODO: `real` only justified for TR
    if !isapprox(χ, χ′, atol=Crystalline.DEFAULT_ATOL)
        error(lazy"expected characters to be convertible to real integers; got $χ")
    end
    χ = χ′

    isₚ, isᵤ = pinned_and_unpinned_operation_indexes(ops)
    nₚ = pinned_symval²ᵀ.(@view ops[isₚ]) # pinned symmetry eigenvalues of 2T

    Nᵤ = length(isᵤ)
    if Nᵤ == 0 # no unpinned symmetry operations
        return Int[], Matrix{Int}(undef, 0, 0), Int[]
    end
    
    χₚ = χ[isₚ, :]
    χᵤ = χ[isᵤ, :]

    χₚᵍ, maxλ_χₚ = generalized_inv(χₚ)
    A = Int.(hcat(χᵤ*(I - χₚᵍ*χₚ/maxλ_χₚ), I))
    Aᵍ, maxλ_A = generalized_inv(A)
    rhs = -Int.((χᵤ*χₚᵍ*nₚ)./maxλ_χₚ)

    fixed = Int.((Aᵍ*rhs)/maxλ_A)
    free  = I - Int.((Aᵍ*A)./maxλ_A)

    opsᵤ = ops[isᵤ]
    fixedᵤ = fixed[end-Nᵤ+1:end]
    freeᵤ = free[end-Nᵤ+1:end, :]
    freeᵤ = foldl(hcat,
                  filter(!iszero, collect(eachcol(freeᵤ))),
                  init=Matrix{Int}(undef, Nᵤ, 0))

    # simplify/reduce freeᵤ: many rows are trivially identical because they associate with
    # operations of the same order - but this challenges lattice reduction techniques, so
    # we manually construct a "operation-order-reduced" (row-wise) matrix, lattice reduce
    # that, and then reconstruct the full `freeᵤ` afterwards
    freeᵤ′ = Matrix{Int}(undef, 0, size(freeᵤ,2))
    treated_rotation_orders = Int[]
    row_positions = Vector{Int}(undef, size(freeᵤ,1))
    for (i, row) in enumerate(eachrow(freeᵤ))
        rot = Crystalline.rotation_order(opsᵤ[i])
        idx = findfirst(==(rot), treated_rotation_orders)
        if idx === nothing
            push!(treated_rotation_orders, rot)
            row_positions[i] = length(treated_rotation_orders)
            freeᵤ′ = vcat(freeᵤ′, row')
        else
            row_positions[i] = idx
            @assert freeᵤ′[idx,:] == row
        end
    end
    @assert freeᵤ == freeᵤ′[row_positions, :]

    if lattice_reduce
        freeᵤ′′ = lattice_reduction(freeᵤ′)
        freeᵤ = freeᵤ′′[row_positions, :] # reconstruct "full" matrix
    end

    return fixedᵤ, freeᵤ, isᵤ
end

# lattice reduce `A` via Seysen's algorithm and canonicalize the output a bit (e.g., first
# row has a nonzero, positive entry in first column))
function lattice_reduction(A)
    isempty(A) && return A
    A′ = seysen(A)[1]
    # to make printing "pretty" and reasonably consistent across cases, we do a few
    # things to normalize `A′` a bit
    if first(A′) == 0
        _, idx = findmax(abs, @view A′[1,:])
        newcols = collect(1:size(A′,2))
        newcols[idx] = 1
        newcols[1] = idx
        A′ = A′[:,newcols]
    end
    if first(A′) < 0
        A′ .*= -1
    end
    return A′
end
# ---------------------------------------------------------------------------------------- #

pinned_symval²ᵀ(op) = round(Int, get_symval²ᵀ(op)) # always integer; enforce in type-domain
function pinned_and_unpinned_operation_indexes(ops::AbstractVector{SymOperation{3}})
    isₚ = Int[] # pinned operator indexes
    isᵤ = Int[] # unpinned operator indexes
    for (i, op) in enumerate(ops)
        rot = Crystalline.rotation_order(op)
        if rot ∈ (-1, -3, -4, -6)
            push!(isᵤ, i)
        else
            push!(isₚ, i)
        end
    end
    return isₚ, isᵤ
end

# ---------------------------------------------------------------------------------------- #

function find_pg_supergroup(pg::PointGroup{D}) where D
    spg = pg # super point group
    sN = Crystalline.order(pg)
    pgiuc = label(pg)
    for pgiuc′ in Iterators.reverse(PG_IUCs[D])
        pgiuc′ == pgiuc && continue
        N′ = Crystalline.PG_ORDERs[pgiuc′]
        # for spg to be a subgroup of pg′, sN must divide N′, cf. Lagrange's theorem
        mod(N′, sN) != 0 && continue 
        pg′ = pointgroup(pgiuc′, Val(D))
        if issubgroup(pg′, pg)
            spg = pg′
            sN = N′
        end
    end
    return spg
end

# ---------------------------------------------------------------------------------------- #

# Computes the generalized inverse `Xᵍ` of `X`, computed from the Smith normal form, 
# multiplied by the maximum absolute value of the Smith norm form singular values, `maxλ`.
# Returns `Xᵍ*maxλ`, which is necessarily an integer matrix, and `maxλ`.
function generalized_inv(X::AbstractMatrix{<:Integer})
    F = smith(X)
    Λg = zeros(Float64, size(diagm(F))[2], size(diagm(F))[1])
    for (n, λₙ) in enumerate(F.SNF)
        Λg[n,n] = iszero(λₙ) ? λₙ : inv(λₙ)
    end
    maxλ = maximum(F.SNF; init=1)
    Xᵍ = round.(Int, F.Tinv*Λg*F.Sinv*maxλ) # generalized inverse × maxλ; necessarily in ℤ

    return Xᵍ, maxλ
end

# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #

# `physical_zero_frequency_gamma_irreps` is the original solution to the problem - but it is
# not well-suited for general-purpose use, since it wants to find the parent point group
# which is both computationally costly & finicky (many corner cases for coordinate system
# differences). But we actually can solve it much simpler by exploiting the insights the
# results from `physical_zero_frequency_gamma_irreps` give, and in particular, the footnote
# in our 2022 PRX paper on O(3) constraints. Working "backwards" from Eq. (S42), we can
# determine a lattice basis for the "free part" of the symmetry eigenvalues, as constrained
# by "integer-ness" of an O(3) irrep expansion. This basis is (see also my summarization in 
# https://math.stackexchange.com/q/5074910/335146):
#   x(-1) = 8∑ₗ cₗ(l+1)
#   x(-3) = 2∑ₗ c₃ₗ - c₃ₗ₊₁
#   x(-4) = 4∑ₗ c₄ₗ - c₄ₗ₊₂
#   x(-6) = 6∑ₗ c₆ₗ + c₆ₗ₊₁ - c₆ₗ₊₃ - c₆ₗ₊₄
# with cₗ ∈ ℤ, and l = 0,1,…,∞. Despite being spanned in an infinite-dimensional space by
# cₗ, the intrinsic dimension of this lattice is 4. The first four "lattice vectors", vᵢ,
# obtained by setting c = [1, 0, 0, 0, 0…], c = [0, 1, 0, 0, 0…], c = [0, 0, 1, 0, 0…], and 
# c = [0, 0, 0, 1, 0…], are linearly independent and span the entire lattice. Thus, the
function _xfree_basis_elements(l)
    v = zeros(Int, 4)
    v[1] = 8 * (l+1)                                       # x(-1)
    v[2] = 2 * (iszero(mod(l,3)) - iszero(mod(l-1,3)))     # x(-3)
    v[3] = 4 * (iszero(mod(l,4)) - iszero(mod(l-2,4)))     # x(-4)
    v[4] = 6 * (iszero(mod(l,6)) + iszero(mod(l-1,6))      # x(-6)
                - iszero(mod(l-3,6)) - iszero(mod(l-4,6)))
    return v
end
const XFREE_BASIS = lattice_reduction(
                        hcat(_xfree_basis_elements(0), _xfree_basis_elements(1),
                             _xfree_basis_elements(2), _xfree_basis_elements(3)))

"""
    physical_zero_frequency_gamma_irreps_O3(
        lgirs :: Collection{LGIrrep{3}}
    )                                       --> Tuple{Vector{Int}, Matrix{Int}}

Equivalent to `physical_zero_frequency_gamma_irreps`, but applies the "super" constraints
from O(3) directly, not merely from the maximal crystallographic super point group.
This implementation is more much more cost-effective, however, and is also more robust to
setting variations. The returned lattice is always reduced via Seysen's algorithm.
Similarly, this implementation always uses the equivalent of `force_fixed = true`.
"""
function physical_zero_frequency_gamma_irreps_O3(lgirs :: Collection{LGIrrep{3}})
    lg = group(lgirs)
    klabel(lg) == "Γ" || error(lazy"input `lgirs` must be Γ-irreps")

    # identify the kinds of roto-inversions that are in the little group
    rotoinv_orders = unique(Crystalline.rotation_order.(lg))
    filter!(∈((-1, -3, -4, -6)), rotoinv_orders) # retain only (non-mirror) roto-inversions
    sort!(rotoinv_orders; rev=true) # sort in (-1, -3, -4, -6) order)
    
    # determine a (reduced) basis for the "free part" of the symmetry eigenvalues, but
    # restricting the lattice to just the symmetry eigenvalues that are actually in `lg`
    include_rows = [r for (o, r) in zip((-1,-3,-4,-6), (1,2,3,4)) if o ∈ rotoinv_orders]
    xfree_basis = XFREE_BASIS[include_rows, :]
    F = smith(xfree_basis)
    d = count(≠(0), F.SNF)
    Λ = Diagonal(@view F.SNF[1:d])
    xfree_basis_projected = F.S*Λ # projected basis for free part of symmetry eigenvalues
    xfree_basis_reduced = lattice_reduction(xfree_basis_projected) # reduced & projected
    @assert size(xfree_basis_reduced) == (length(rotoinv_orders), length(rotoinv_orders))

    # express the symmetry eigenvalues as a lattice with fixed component `xfixed` and a
    # lattice basis term `xfree`, i.e., as `x = xfixed + xfree * [l, p, q, …]`
    xfixed = zeros(Int, length(lg))
    xfree = zeros(Int, length(lg), length(rotoinv_orders))
    for (row, op) in enumerate(lg)
        xfixed[row] = pinned_symval²ᵀ(op)
        length(rotoinv_orders) == 0 && continue
        
        # set `xfree[row, :]` according to the corresponding row of `xfree_basis_reduced`;
        # we just need to determine what this row is:
        row′ = findfirst(==(Crystalline.rotation_order(op)), rotoinv_orders)
        isnothing(row′) && continue # not a roto-inversion operation; everything pinned
        xfree[row, :] .= @view xfree_basis_reduced[row′, :]
    end

    # now we have the symmetry eigenvalues expressed as a lattice, and we just need to
    # convert this to a lattice in the irrep multiplicities of `lgirs`
    ct = characters(lgirs)
    nfixed = ct\xfixed
    nfixed_i = round.(Int, real.(nfixed))
    if !isapprox(nfixed, nfixed_i, atol=1e-10)
        error(lazy"expected real, integer solutions for `nfixed`, got $nfixed")
    end
    nfree = ct\xfree
    nfree_i = round.(Int, real.(nfree))
    if !isapprox(nfree, nfree_i, atol=1e-10)
        error(lazy"expected real, integer solutions for `nfree`, got $nfree")
    end

    # we can return immediately; there's no permutations to consider, since we never went
    # away from the original `lgirs` (i.e., we didn't involve another concrete point group)
    return nfixed_i, nfree_i
end
