using Crystalline
using Crystalline: smith, classification, prettyprint_symmetryvector
using PhotonicBandConnectivity
using SymmetryBases
using PrettyTables
using Test

const PBC = PhotonicBandConnectivity

include("text_utils.jl"); using Main.TextUtils

#------------------------------------------------------------------------------------------
# SETUP

sgnums       = 1:230
timereversal = false
Nsolutions   = 3 # how many μᵀ solutions we want
μᵗstart      = 3 # start target filling
verbose      = false

for sgnum in sgnums
    print("SG $(sgnum): ")

    # -------------------------------------------------------------------------------------
    # PREP-WORK (IRREPS, HILBERT BASES, PINNED EIGENVALUES)

    lgirsd = lgirreps(sgnum, Val(3))
    timereversal && (foreach(zip(keys(lgirsd), values(lgirsd))) do (klab, lgirs)
                        lgirsd[klab] = realify(lgirs)
                     end)
    lgirs_Γ = lgirsd["Γ"]

    # prep-work to get Hilbert bases etc
    brs  = bandreps(sgnum, spinful=false, timereversal=timereversal)
    B    = stack(brs) # matrix with columns of EBRs.
    isℤ₁ = classification(brs) == "Z₁"
    if isℤ₁
        println(" ... skipping; trivial symmetry indicator group ℤ₁\n") 
        continue
    end
    println()

    F    = smith(B)   # Smith normal decomposition of B
    Nⁱʳʳ = size(B, 1) # number of irreps plus 1 (filling)
    sb   = compatibility_basis(F, brs)
    Γidxs    = PBC.get_Γidxs(lgirs_Γ, sb)   
    notΓidxs = [idx for idx in 1:Nⁱʳʳ if idx ∉ Γidxs]

    ms¹ᴸ = PBC.find_representation¹ᴸ(lgirs_Γ)
    ms²ᵀ = PBC.find_representation²ᵀ(lgirs_Γ)
    ms   = PBC.find_representation²ᵀ⁺¹ᴸ(lgirs_Γ)

    ntidxs¹ᴸ  = PBC.find_symmetry_constrained_bases(sb, ms¹ᴸ, Γidxs)
    μsᴴ       = PBC.fillings(sb)
    _, pick¹ᴸ = findmin(μsᴴ[ntidxs¹ᴸ])
    idx¹ᴸ     = ntidxs¹ᴸ[pick¹ᴸ]

    #--------------------------------------------------------------------------------------
    # FINDING VALID SOLUTIONS

    μᵗ               = μᵗstart-1
    Nsolutions_found = 0
    Γidx_in_sb       = findfirst(==("Γ"), sb.klabs)
    
    while true
        flush(stdout)

        μᵗ += 1
        cⁱs, μᵀ = check_target_filling_regular1L(μᵗ, ms¹ᴸ, ms, μsᴴ, sb, idx¹ᴸ, Γidxs, 
                                                 notΓidxs; verbose=false)
        
        isempty(cⁱs) && continue # go to next μᵗ if no solutions found

        # Extract the "physical parts" of the symmetry vector (i.e. sans singular Γ-irreps)
        nᵀ⁺ᴸs = PBC.sum_symbases.(Ref(sb), cⁱs)
        nᴸ    = sb[idx¹ᴸ]
        nᵀs   = nᵀ⁺ᴸs .- Ref(nᴸ)
        for nᵀ in nᵀs
            # tricky indexing: this achieves the ordering-adjusted subtraction into the `sb`
            # basis (note the broadcasted assignment, which is what makes this work)
            nᵀ[Γidxs] .-= ms²ᵀ
        end
        unique!(nᵀs) # remove equivalent solutions (due to expansion non-uniqueness)
        
        # construct human-readable symmetry vectors
        ios = [IOBuffer() for _ in eachindex(nᵀs)]
        prettyprint_symmetryvector.(ios, nᵀs, Ref(sb.irlabs))
        nᵀs_str = String.(take!.(ios))

        # insert (▪)²ᵀ in Γ-irrep spot
        nᵀs_str .= add_content_to_symvec_str_at_kidx.(nᵀs_str, Γidx_in_sb, Ref("(▪)²ᵀ")) # from /test/text_utils.jl

        # write results to screen
        minN   = Nsolutions_found + 1
        plural = length(nᵀs) > 1

        println(
            "   ", minN, ordinal_indicator(minN)," minimal solution: ",
            "μᵀ = ", μᵀ, " (", length(nᵀs), " symmetry vector", plural ? "s" : "", ")")
        flush(stdout)

        # find "Z₂" factor-type topology of each solution
        topos = Vector{TopologyKind}(undef, length(nᵀs))
        filling_enforced = true
        for (idx, nᵀ) in enumerate(nᵀs)
            topos[idx] = calc_topology_singular(nᵀ, nᴸ, ms²ᵀ, Γidxs, F)
            topos[idx] == NONTRIVIAL || (filling_enforced = false; break)
        end

        if filling_enforced
            println("      ⇒ Filling-enforced topology!")
            pretty_table([nᵀs_str topos], ["nᵀ", "topology"]; # contents & header row
                            alignment = :l, crop = :none, tf = tf_unicode, 
                            vlines = :none, hlines = [:begin, 1, :end])
        end

        # check if we are done
        Nsolutions_found += 1
        Nsolutions_found ≥ Nsolutions && break
    end # while
    println()
end # for sgnum