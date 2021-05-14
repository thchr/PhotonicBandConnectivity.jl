using Pkg
dirname(Pkg.project().path) == (@__DIR__) || (Pkg.activate(@__DIR__); cd(@__DIR__))

using Crystalline
using Crystalline: prettyprint_symmetryvector
using SymmetryBases
using PhotonicBandConnectivity
using PhotonicBandConnectivity: minimal_expansion_of_zero_freq²ᵀ⁺¹ᴸ_bands
using ProgressMeter
using JLD2

# find all the minimum fillings below the fundamental gap in systems where the longitudinal
# bands are real, such as in phononic or magnonic systems
sgnums   = 1:230
has_tr   = false
savedata = false
show_has_fragile_phases = false

μ²ᵀ⁺¹ᴸs = Vector{Int}(undef, length(sgnums)) # minimal fillings below fundamental gap
μs      = Vector{Int}(undef, length(sgnums))
topos   = Vector{Vector{TopologyKind}}(undef, length(sgnums))
for (sgidx, sgnum) in enumerate(sgnums)
    print("SG ", sgnum, ":\t")
    
    cⁱs, μ²ᵀ⁺¹ᴸs[sgidx], sb = minimal_expansion_of_zero_freq²ᵀ⁺¹ᴸ_bands(sgnum, has_tr)
    μs[sgidx] = minimum(fillings(sb)) # "default" fillings
    print("μ²ᵀ⁺¹ᴸ = ", μ²ᵀ⁺¹ᴸs[sgidx], " (μ = ", μs[sgidx], ")")

    # calculate unique symmetry vectors
    n²ᵀ⁺¹ᴸs = unique!(sort(sum_symbases.(Ref(sb), cⁱs)))
    print(" w/ ", length(n²ᵀ⁺¹ᴸs), " solution", length(n²ᵀ⁺¹ᴸs)≠1 ? "s" : "")

    # get EBRs
    BRS = bandreps(sgnum, 3, timereversal=has_tr)
    B   = matrix(BRS, true)
    F   = Crystalline.smith(B)

    # print symmetry-identifiable class of space group
    print("\t [", classification(BRS))
    if show_has_fragile_phases 
        if !has_tr && sgnum ∈ (83,174,175,176)
            # the non-topological basis for tr-broken SGs 83, 174, 175, 176 is extremely
            # slow (does not terminate in +24h); so we just don't bother to print the
            # general information about fragility here
            print(" | skipped fragility check (quagmire)]")
        else
            nontopo_sb = nontopological_basis(F, BRS)
            trivial_idxs, fragile_idxs = split_fragiletrivial(nontopo_sb, B)
            !isempty(fragile_idxs) && print("+fragile")
        end
    end
    println("]"); flush(stdout)

    # calculate trivial/nontrivial/fragile topology classification of all solutions
    topos[sgidx] = @showprogress map(n -> calc_detailed_topology(n, B, F), n²ᵀ⁺¹ᴸs)

    # print if filling-enforced topology detected
    if all(==(NONTRIVIAL), topos[sgidx])
        println("\t⇒ Filling-enforced stable topology!")
    elseif all(==(FRAGILE), topos[sgidx])
        println("\t⇒ Filling-enforced fragile topology!")
    elseif all(≠(TRIVIAL), topos[sgidx])
        println("\t⇒ Filling-enforced mixed stable & fragile topology")
    end
end

# save data, if requested
if savedata && sgnums == 1:230 && has_tr
    @save "phonon-connectivity-result-with-tr.jld2" sgnums has_tr μ²ᵀ⁺¹ᴸs μs topos
elseif savedata && sgnums == 1:230 && !has_tr
    @save "phonon-connectivity-result-without-tr.jld2" sgnums has_tr μ²ᵀ⁺¹ᴸs μs topos
end