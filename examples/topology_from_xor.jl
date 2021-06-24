using PhotonicBandConnectivity
using SymmetryBases
using Test
using Crystalline
using Crystalline: matrix

#= 
This is a test script to check the topology of the 2T solutions in the cases where the 2T 
solutions is irregular, such that it requires inclusion of a 1L mode. The approach is to 
contrast the topology of the 1L mode choice/pick with the topology of the 2T+1L solution.
See src/topology_as_2T1L_vs_1L_difference.jl for additional details.
=#

sgnums = 1:230
has_tr = false # true

# ω=0 solutions
data   = minimal_expansion_of_zero_freq_bands.(sgnums, timereversal=has_tr);
cⁱss   = getindex.(data, 1) # coefficients of expansions
νᵀs    = getindex.(data, 2) # fillings for tranverse branch
sbs    = getindex.(data, 3) # symmetry bases
idx¹ᴸs = getindex.(data, 4) # index for chosen 1L branch

# Band representations (to check trivial vs. nontrivial topology)
BRSs = bandreps.(sgnums, timereversal=has_tr)

# check whether 1L pick is trivial
for (sgidx, sgnum) in enumerate(sgnums)
    idx¹ᴸ = idx¹ᴸs[sgidx]
    
    # TODO: Treat the simpler `idx¹ᴸ !== nothing` case as well within this script?
    if idx¹ᴸ !== nothing && classification(BRSs[sgidx]) ≠ "Z₁"
        sb, BRS = sbs[sgidx], BRSs[sgidx]
        B  = matrix(BRS)
        F  = smith(B)
        nᴸ = sb[idx¹ᴸ]

        # -------------------------------------------------------------------------------- #
        # check topology of all T+L solutions using `topology_from_2T1L_xor_1L` from
        # src/topology_as_2T1L_vs_1L_difference.jl (xor difference of T+L and L topology)

        # find unique symmetry vectors
        ns = unique!(sort(sum_symbases.(Ref(sb), cⁱss[sgidx])))

        println("SG ", sgnum, ": νᵀ = ", νᵀs[sgidx], 
                " (", length(ns), " ω=0 solution", length(ns) > 1 ? "s)" : ")")

        # find "Z₂" factor-type topology of each solution (this step can be a little slow 
        # for some of the SGs with many solutions/large M, because the optimization step 
        # is a bit slow)
        topos = topology_from_2T1L_xor_1L.(ns, Ref(nᴸ), Ref(F))

        # get aggregated stats 
        trivial_countᵀ    = count(==(TRIVIAL),    topos)
        nontrivial_countᵀ = count(==(NONTRIVIAL), topos)

        # print summary of transverse solutions' topology
        tdigs, ntdigs = ndigits(trivial_countᵀ), ndigits(nontrivial_countᵀ)
        maxdigs = max(tdigs, ntdigs)
        println("   # trivial:    ", " "^(maxdigs-tdigs),  trivial_countᵀ)
        println("   # nontrivial: ", " "^(maxdigs-ntdigs), nontrivial_countᵀ)
        trivial_countᵀ == 0 && println("   ─── all-nontrivial solution ───")
        #println()
        # NOTE: SGs 13, 48, 50, 68, 86 are interesting (all transverse ω=0 solutions nontrivial)

        # --- check topology via "index-difference" ---
        indices′_and_Λ = indicators_singular.(ns, Ref(nᴸ), Ref(F))
        indices′ = first.(indices′_and_Λ)
        Λ = first(last.(indices′_and_Λ))

        topos′ = ifelse.(iszero.(indices′), TRIVIAL, NONTRIVIAL)

        # get aggregated stats 
        trivial_countᵀ′    = count(==(TRIVIAL),    topos′)
        nontrivial_countᵀ′ = count(==(NONTRIVIAL), topos′)

        # print summary of transverse solutions' topology
        tdigs′, ntdigs′ = ndigits(trivial_countᵀ′), ndigits(nontrivial_countᵀ′)
        maxdigs′ = max(tdigs′, ntdigs′)
        println("   # trivial:    ", " "^(maxdigs′-tdigs′),  trivial_countᵀ′)
        println("   # nontrivial: ", " "^(maxdigs′-ntdigs′), nontrivial_countᵀ′)
        trivial_countᵀ′ == 0 && println("   ─── all-nontrivial solution ───")
        

        topos ≠ topos′ && @info "Difference! Λ = $Λ"
        topos ≠ topos′ && @info "Indices: " indices′
        topos ≠ topos′ && trivial_countᵀ′ == 0 && @warn "Filling-enforced topology!"
        println()
    end
end

# TODO: we should check that we get the same T Z₂ index by computing whether or not we can 
#       expand the non-Γ parts of the T solutions in a strictly positive coefficient EBR 
#       expansion