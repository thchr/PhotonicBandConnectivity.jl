module PhotonicBandConnectivity

using Crystalline
using Crystalline: rotation, AbstractSymmetryVector
using Crystalline: bandreps, Collection, BandRep, LGIrrep, SiteIrrep
using SymmetryBases
using SymmetryBases: PyNormaliz
using PythonCall
using LinearAlgebra
using DocStringExtensions

export minimal_expansion_of_zero_freq_bands,
    sum_symbases, sum_symbases!,
    check_target_filling_regular1L,
    calc_topology_singular,
    indicators_singular,
    topology_from_2T1L_xor_1L, # deprecated; remove eventually
    is_transverse_bandstruct,
    transverse_symmetry_vectors,
    transverse_vrep,
    is_vrep

# ---------------------------------------------------------------------------------------- #

# alias for a spinless band representation; spinful variants are irrelevant for photons
const SpinlessBandRep{D} = BandRep{D, LGIrrep{D}, SiteIrrep{D}}

# ---------------------------------------------------------------------------------------- #

include("planewave_symvals.jl")
include("irreps_and_representations.jl")
include("symbasis_utils.jl")
include("constrained_expansions.jl")
#include("src/legacy_constrained_expansions.jl")
include("topology_as_2T1L_vs_1L_difference.jl")
include("is_transverse_bandstruct.jl")

include("transverse_symmetry_vector_solutions.jl")
include("transverse_vrep.jl")

# ---------------------------------------------------------------------------------------- #
# TODO: Use the following to make the output of `minimal_expansion_of_zero_freq_bands` more
#       structured, user-friendly, and sane

#=
function minimal_transverse_vectors(
    sgnum::Integer, Dᵛ::Val{D}=Val(3);
    timereversal::Bool=true
    ) where D
    
    D ≠ 3 && error("only D = 3 is supported; D = 2 can be treated as regular")

    lgirsd = lgirreps(sgnum, Dᵛ)
    timereversal && realify!(lgirsd)
    lgirsΓ = lgirsd["Γ"]
    vrep = transverse_vrep(lgirsΓ)
    push!(parent(lgirsΓ), vrep) # append fake vrep to `lgirsΓ`
    vrep_lab = label(vrep)

    cⁱs, _, sb, idxᴸ = minimal_expansion_of_zero_freq_bands(sgnum; timereversal)
    Nⁱʳ = length(first(sb))
    _ns = map(cⁱs) do cⁱ
        _n = zeros(Int, Nⁱʳ)
        sum_symbases!(_n, sb, cⁱ)
        !isnothing(idxᴸ) && (_n .-= sb[idxᴸ])
        _n
    end
    lastΓidx = something(findlast(irlab->klabel(irlab)=="Γ", sb.irlabs))     

    # forcibly add the virtual rep to the symmetry vectors (we subtract double-count next)
    foreach(_ns) do n
        insert!(n, lastΓidx+1, 1)
    end
    irlabs = insert!(copy(sb.irlabs), lastΓidx+1, vrep_lab)
    ns = SymmetryVector.(_ns, Ref(irlabs), Ref(lgirsd))

    # at this point, we're double-counting the transverse virtual rep; remove it from 
    # the regular irrep multiplicities
    lgirsv_in_ns = first(ns).lgirsv
    Γidx = something(findfirst(lgirs->klabel(first(lgirs))=="Γ", lgirsv_in_ns))
    lgirsΓ_in_ns = first(ns).lgirsv[Γidx] # may have different sorting from `lgirsd[Γ]` (but vrep is still last)
    cs = find_representation²ᵀ(@view lgirsΓ_in_ns[1:end-1])
    for n in ns
        mults = multiplicities(n)[Γidx]
        mults[1:end-1] .-= cs
        @assert all(≥(0), mults)
    end

    return ns, lgirsΓ[end], lgirsd
end
=#

# ---------------------------------------------------------------------------------------- #

function minimal_expansion_of_zero_freq_bands(
            sgnum::Integer; 
            timereversal::Bool=true, verbose::Bool=false, shuffle_1Lpick::Bool=false,
            allpaths::Bool=false, safetychecks::Bool=true
            )

    # Irreps at Γ, irrep-multiplicities of ω=0 2T bands, and symmetry operations
    lgirs = lgirreps(sgnum, Val(3))["Γ"]
    timereversal && (lgirs = realify(lgirs))
    lg = group(first(lgirs))
    rotvals = map(op->Crystalline.rotation_order(op), lg)

    # 2T irreps; check if "simple treatment"/fast-path is applicable
    ms²ᵀ = find_representation²ᵀ(lgirs)
    has_nonmirror_improper = any(∈((-1, -3, -4, -6)), rotvals)
    is_regular²ᵀ = all(≥(0), ms²ᵀ)
    if !has_nonmirror_improper && is_regular²ᵀ 
        # Scenario (easy case): All symvals known & regular 2T irrep

        # If there are no non-mirror improper rotations, we can directly infer the irrep of
        # the 2T branches. If that irrep is regular (i.e. has no negative coefficients), we
        # don't need to invoke 1L at all, and can solve for just 2T alone.
        r = find_minimum_bandreps_regular2T(sgnum, lgirs, timereversal, ms²ᵀ; 
                                            safetychecks=safetychecks, verbose=verbose,
                                            allpaths=allpaths)
        return r..., nothing
        
    else 
        # Two possible scenarios (treat with same approach):
        #   - All symvals known & irregular 2T and regular 1L irreps
        #   - Not all symvals known; multiple irrep options
        ms¹ᴸ = find_representation¹ᴸ(lgirs)
        ms   = ms²ᵀ .+ ms¹ᴸ
        @assert all(ms .== ms²ᵀ .+ ms¹ᴸ)
        @assert all(≥(0), ms¹ᴸ)                        # check: 1L irrep regular (Γ₁)
        @assert ms == find_representation²ᵀ⁺¹ᴸ(lgirs)  # →: [2T+1L] = 2T+1L

        return find_minimum_bandreps_regular1L(sgnum, lgirs, timereversal, ms¹ᴸ, ms;
                                               verbose=verbose, 
                                               shuffle_1Lpick=shuffle_1Lpick,
                                               allpaths=allpaths)
    end
    # TODO: The returned coefficients cⁱs do not necessarily each describe different 
    #       symmetry vectors n, since the Hilbert basis is not linearly independent.
    #       We should consider whether it would be better to only return expansions that 
    #       result in unique symmetry vectors, i.e. in general a subset of cⁱs
end

function find_minimum_bandreps_regular2T(sgnum, lgirs, timereversal, ms²ᵀ; 
                                         verbose::Bool=false, safetychecks::Bool=false,
                                         allpaths::Bool=false)
    verbose && println("SG ", sgnum)

    sb, Γidxs = compatibility_basis_and_Γidxs(lgirs; timereversal, allpaths)
    νsᴴ = fillings(sb)
    νᴴₘₐₓ = maximum(νsᴴ)

    # We seek an expansion with coefficients cᵢ≥0 such that
    #   P(Γ) ∑ᵢ cᵢ 𝐧ᴴᵢ ≥ 𝐦(Γ)
    # where P(Γ) projects out the Γ-irreps from the Hilbert bases 𝐧ᴴᵢ. In code, letting
    # `nsᴴ = stack(sb)`, this means we seek a solution with `nsᴴ[Γidxs,:]*c ≥ ms`. 
    # Finally, we impose a filling
    # constraint, such that the overall number of bands is at most ν. In code, this requires
    # that `nsᴴ[end,:]*c == ν`. Moreover, all 𝐧ᴴᵢ that does not have at least one nonzero
    # entry matching `ms` will not help us in fulfilling these constraints in a nontrivial
    # way, so we can ignore those (would correspond to just stacking on some bands).
    # Finally, we can restrict the sum to contain at most two 𝐧ᴴᵢ (same or different): if we
    # have more, then at least one of them isn't actually needed to fulfil ``𝐧(Γ) ≥ 𝐦(Γ)``,
    # and can then be considered a trivial stacking.

    # the "nontrivial" parts of `nᴴ` must have at least one positive element for the same 
    # irrep as a nonzero index of `ms`; we can ignore all the others
    ntidxs²ᵀ = find_symmetry_constrained_bases(sb, ms²ᵀ, Γidxs)

    maxdepth = 2
    for ν²ᵀᵗ in 2:2νᴴₘₐₓ # target filling for 2T branches (≥2)
        cⁱs = filling_symmetry_constrained_expansions(ν²ᵀᵗ, ms²ᵀ, νsᴴ, sb, Γidxs;
                                    ntidxs = ntidxs²ᵀ,   # include only "nontrivial" bases
                                    maxdepth = maxdepth) # limit to two basis terms

        if !isempty(cⁱs)
            verbose && println("   ⇒ νᵀ = ", ν²ᵀᵗ, ": ", length(cⁱs), " valid expansions")            
            safetychecks && safetycheck²ᵀ(cⁱs, ν²ᵀᵗ, ms²ᵀ, νsᴴ, sb, Γidxs)
            
            return cⁱs, ν²ᵀᵗ, sb
        end
    end
    throw("Found no valid expansions consistent with constraints")
end

function find_minimum_bandreps_regular1L(sgnum, lgirs, timereversal, ms¹ᴸ, ms;
                                         verbose::Bool=false, shuffle_1Lpick::Bool=false,
                                         allpaths::Bool=false)
    verbose && print("SG ", sgnum)

    sb, Γidxs = compatibility_basis_and_Γidxs(lgirs; timereversal, allpaths)
    Nⁱʳʳ = length(first(sb))
    notΓidxs = [idx for idx in 1:Nⁱʳʳ if idx ∉ Γidxs]
    νsᴴ = fillings(sb)
    νᴴₘᵢₙ, νᴴₘₐₓ = extrema(νsᴴ)
    verbose && println(" (νᴴₘᵢₙ = ", νᴴₘᵢₙ, ")")
    
    # Here, the irrep of 1L is regular (Γ₁) and the irrep of 2T is irregular (i.e. has 
    # negative coefficients). As a result, it is impossible to expand 2T's irrep in the
    # Hilbert basis since it has strictly positive elements and coefficients. We can still
    # can try to find an expansion for 2T+1L simultaneously.
    ntidxs¹ᴸ  = find_symmetry_constrained_bases(sb, ms¹ᴸ, Γidxs)
    _, pick¹ᴸ = findmin(νsᴴ[ntidxs¹ᴸ])
    if shuffle_1Lpick 
        if length(ntidxs¹ᴸ) > pick¹ᴸ
            pick¹ᴸ += 1
        else
            @warn "Unable to shuffle 1L pick"
        end
    end
    idx¹ᴸ = ntidxs¹ᴸ[pick¹ᴸ] # I've tested that resulting expansions for 2T (=[1L+2T]-[1L]) are invariant wrt. to this choice
    νᴸ = νsᴴ[idx¹ᴸ]
        
    # find _all_ feasible solutions to ms constraints for fixed and minimal νᵗ; we use a 
    # recursive looping construct to find candidate expansions
    max_patience_νᵗ = max(4*νᴴₘₐₓ, 12)
    for νᵗ in (2+νᴸ):max_patience_νᵗ # target filling (≥2+νᴸ) (function returns from loop)
        verbose && print("   … νᵗ = ", νᵗ, ": ")

        cⁱs, νᵀ = check_target_filling_regular1L(νᵗ, ms¹ᴸ, ms, νsᴴ, sb, idx¹ᴸ, Γidxs, notΓidxs;
                                                 verbose=verbose)

        if !isempty(cⁱs)
            verbose && println("   ⇒ νᵀ = ", νᵀ, ": ", length(cⁱs), " valid expansions")
            return cⁱs, νᵀ, sb, idx¹ᴸ
        end
    end
    
    throw("Found no valid expansions consistent with constraints")
end

# Note that this can be used generically for both 2T and 1L cases
"""
    $(SIGNATURES)

For a total (T+L) target filling `νᵗ`, subject to Γ ``ω=0`` constraints `ms¹ᴸ` and `ms`, 
determine whether a solution to the symmetry constraints exist in the space group with
compatibility (Hilbert) basis `sb::SymBasis` [with associated fillings `νsᴴ=fillings(sb)`].
If a solution exists, returns the associated solutions as vectors `cⁱs` indexing into `sb` 
as well as the associated transverse filling `νᵀ`. The latter is inferred from the total
filling `νᵗ` and the filling of the longitudinal (1L) mode choice, `νᴸ` (which is indicated
by the 1L index `idx¹ᴸ`, which gives the longitudinal symmetry vector `nᴸ = sb[idx¹ᴸ]`).
"""
function check_target_filling_regular1L(νᵗ, ms¹ᴸ, ms, νsᴴ, sb::SymBasis, idx¹ᴸ, 
            Γidxs, notΓidxs; verbose::Bool=false)
    # Find the solutions to c₁νᴴ₁ + c₂νᴴ₂ + ... = νᵗ subject to the 2T+1L ms constraint
    # The below basically uses recursion to do a nested set `max_terms` loops which
    # solves linear Diophantine equation and checks symmetry constraints as well; the
    # maximum number of included bases in a valid expansion is div(νᵗ, νᴴₘᵢₙ, RoundDown)
    cⁱs_candidates = filling_symmetry_constrained_expansions(νᵗ, ms, νsᴴ, sb, Γidxs)
    verbose && println(length(cⁱs_candidates), " expansion candidates")

    nᴸ = sb[idx¹ᴸ]  # symmetry vector of 1L pick
    νᴸ = νsᴴ[idx¹ᴸ] # assoc. longitudinal filling
    νᵀ = νᵗ - νᴸ    # assoc. transverse filling
    
    # Proceed to check combinations of nᴸ and n=sum(sb[cⁱ])
    cⁱs = Vector{Int}[]
    n = similar(first(sb)) # candidate solution buffer     
    for cⁱ in cⁱs_candidates # 2T+1L constraints
        sum_symbases!(n, sb, cⁱ) # compute new candidate vector from cⁱ indices
        # test 1: n(∉Γ)-nᴸ(∉Γ) ≥ 0
        if all(≥(0), @views n[notΓidxs] .- nᴸ[notΓidxs])
            # test 2: [n(Γ)-n¹ᴸ⁺²ᵀₚᵢₙ(Γ)] - [nᴸ(Γ)-n¹ᴸₚᵢₙ(Γ)] ≥ 0 
            # note: this is the same as nᵀ(Γ) - n²ᵀₚᵢₙ(Γ) ≥ 0 where nᵀ ≡ n-nᴸ and
            #       n²ᵀₚᵢₙ(Γ) ≡ n¹ᴸ⁺²ᵀₚᵢₙ(Γ)-n¹ᴸₚᵢₙ(Γ). I.e. just the non-singular/physical
            #       Γ irreps of higher-lying bands not connected to ω=0 directly
            if all(≥(0), (n[Γidxs] .- ms) .- (nᴸ[Γidxs] .- ms¹ᴸ)) 
                push!(cⁱs, cⁱ) # found a valid solution; push to storage
            end
        end
    end

    return cⁱs, νᵀ    # if no solutions were valid, `cⁱs` will be empty
end

# ---------------------------------------------------------------------------------------- #
# A simpler version of the above problem is to calculate the minimal connectivities in
# a phononic or magnonic problem, where the 1L modes are real. Then we don't need to
# subsequently factor out the 1L solutions - instead, we can work with the 1L+2L solely.
# The function below implements this, returning the minimal connectivities in bosonic
# systems where the longitudinal modes are real.

function minimal_expansion_of_zero_freq²ᵀ⁺¹ᴸ_bands(sgnum::Integer;
            timereversal::Bool=true, verbose::Bool=false, allpaths::Bool=false, 
            safetychecks::Bool=false)

    # Irreps at Γ, irrep-multiplicities of ω=0 2T bands, and symmetry operations
    lgirs   = lgirreps(sgnum, Val(3))["Γ"]
    timereversal && (lgirs = realify(lgirs))
    lg      = group(first(lgirs))

    # Hilbert bases and "default" connectivities
    sb, Γidxs = compatibility_basis_and_Γidxs(lgirs; timereversal, allpaths)
    νsᴴ       = fillings(sb)
    νᴴₘₐₓ     = maximum(νsᴴ)

    # 2T+1L symmetry constraints
    ms     = find_representation²ᵀ⁺¹ᴸ(lgirs)
    ntidxs = find_symmetry_constrained_bases(sb, ms, Γidxs)

    for νᵗ in 3:3νᴴₘₐₓ # target filling for 2T+1L branches (≥3)
        cⁱs = filling_symmetry_constrained_expansions(νᵗ, ms, νsᴴ, sb, Γidxs; 
                                    ntidxs = ntidxs)  # include only "nontrivial" bases

        if !isempty(cⁱs)
            verbose && println("   ⇒ νᵀ = ", νᵗ, ": ", length(cⁱs), " valid expansions")            
            safetychecks && safetycheck²ᵀ(cⁱs, νᵗ, ms, νsᴴ, sb, Γidxs)
            
            return cⁱs, νᵗ, sb
        end
    end
    throw("Found no valid expansions consistent with constraints")
end

# ---------------------------------------------------------------------------------------- #

include("unpinned_gamma_irrep_solution_space.jl")
export physical_zero_frequency_gamma_irreps_O3
export prettyprint_irrep_solspace

# ---------------------------------------------------------------------------------------- #

end # module