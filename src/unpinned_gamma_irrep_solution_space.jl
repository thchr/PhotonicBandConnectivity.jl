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
            lattice_reduce         :: Bool = true
            )
    
    if any(lgir->reality(lgir) ≠ REAL && !lgir.iscorep, lgirs)
        # non-`realified` irreps; presently not handled
        error("""input `lgirs` must be a set of `realified` or inherently real irreps;
                 time-reversal breaking not yet supported""")                
    end
    lg = group(first(lgirs))
    if klabel(lg) != "Γ" # only Γ-irreps are supported
        error(lazy"input `lgirs`` are associated with the $(label(lg))-point: must be Γ")
    end
    pg = something(find_parent_pointgroup(lg))
    pgirs = realify(pgirreps(label(pg))) # TODO: Generalize to non-`realified` irreps
    pgops = operations(group(first(pgirs)))

    local fixed, free, unpinned_idxmap
    if supergroup_constraints
        spg = find_pg_supergroup(pg)
        spgirs = realify(pgirreps(label(spg)))
        spgops = operations(group(first(spgirs)))
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

    # now we have a representation of the symeigs that we like; convert it to a irrep 
    # multiplicities
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
        freeᵤ′′ = seysen(freeᵤ′)[1] # lattice reduce via Seysen algorithm
        # now, to make printing "pretty" and reasonably consistent across cases, we do a few
        # things to normalize the output a bit
        if first(freeᵤ′′) == 0
            _, idx = findmax(abs.(freeᵤ′′[1,:]))
            newcols = collect(1:size(freeᵤ′′,2))
            newcols[idx] = 1
            newcols[1] = idx
            @assert sort(newcols) == 1:size(freeᵤ′′,2)
            freeᵤ′′ = freeᵤ′′[:,newcols]
        end
        if first(freeᵤ′′) < 0
            freeᵤ′′ .*= -1
        end

        freeᵤ = freeᵤ′′[row_positions, :] # reconstruct "full" matrix
    end

    return fixedᵤ, freeᵤ, isᵤ
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