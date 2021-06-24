using Crystalline
using Crystalline: matrix, smith, classification, prettyprint_symmetryvector
using PhotonicBandConnectivity
using SymmetryBases
using PrettyTables
using Test

const PBC = PhotonicBandConnectivity

includet("text_utils.jl"); using Main.TextUtils

#------------------------------------------------------------------------------------------
# TABLE CONFIG
function table_config(type::String)
    if type == "latex"
        return (backend = :latex, hlines= Int[], 
                table_type = :longtable, 
                longtable_footer = "\\emph{\\footnotesize\\ldots\\ continued on next page}")
    elseif type == "unicode"
        return (crop = :none, tf = tf_unicode, vlines = :none, hlines = [:begin, 1, :end])
    elseif type == "markdown"
        return (crop = :none, tf = tf_unicode)
    end
end


#------------------------------------------------------------------------------------------
# SETUP

sgnum        = 87#1:230#1:230
timereversal = true
Nsolutions   = 1 # how many νᵀ solutions we want
νᵗstart      = 2 # start target filling
outputtype   = "unicode"
io           = stdout
verbose      = true
t₀           = time()

repr²ᵀ = outputtype ∈ ("markdown", "unicode") ? "(▪)²ᵀ" : "(\\smallsquare)\\tT"


println(io, outputtype=="latex" ? "\\subsubsection*{SG $sgnum}\n" : "## SG $(sgnum)\n")
io ≠ stdout && println("SG $(sgnum)")

# -------------------------------------------------------------------------------------
# PREP-WORK (IRREPS, HILBERT BASES, PINNED EIGENVALUES)

lgirsd = get_lgirreps(sgnum, Val(3))
timereversal && (foreach(zip(keys(lgirsd), values(lgirsd))) do (klab, lgirs)
                    lgirsd[klab] = realify(lgirs)
                    end)
lgirs_Γ = lgirsd["Γ"]

# prep-work to get Hilbert bases etc
BRS  = bandreps(sgnum, spinful=false, timereversal=timereversal)
B    = matrix(BRS) # Matrix with columns of EBRs.
F    = smith(B)    # Smith normal decomposition of B
Nⁱʳʳ = size(B, 1)  # number of irreps plus 1 (filling)
isℤ₁ = classification(BRS) == "Z₁"
    
sb       = compatibility_basis(F, BRS)
Γidxs    = PBC.get_Γidxs(lgirs_Γ, sb)   
notΓidxs = [idx for idx in 1:Nⁱʳʳ if idx ∉ Γidxs]
io ≠ stdout && println("   ... found general Hilbert basis (", length(sb), " bases)")

ms¹ᴸ = PBC.find_representation¹ᴸ(lgirs_Γ)
ms²ᵀ = PBC.find_representation²ᵀ(lgirs_Γ)
ms   = PBC.find_representation²ᵀ⁺¹ᴸ(lgirs_Γ)

ntidxs¹ᴸ  = PBC.find_symmetry_constrained_bases(sb, ms¹ᴸ, Γidxs)
νsᴴ       = PBC.fillings(sb)
_, pick¹ᴸ = findmin(νsᴴ[ntidxs¹ᴸ])
idx¹ᴸ     = ntidxs¹ᴸ[pick¹ᴸ]

#--------------------------------------------------------------------------------------
# FINDING VALID SOLUTIONS

νᵗ               = νᵗstart-1
Nsolutions_found = 0
Γidx_in_sb       = findfirst(==("Γ"), sb.klabs)

while true
    global νᵗ, Nsolutions_found
    νᵗ += 1
    cⁱs, νᵀ = check_target_filling_regular1L(νᵗ, ms¹ᴸ, ms, νsᴴ, sb, idx¹ᴸ, Γidxs, 
                                                notΓidxs; verbose=verbose)
    
    isempty(cⁱs) && continue # go to next νᵗ if no solutions found

    # Extract the "physical parts" of the symmetry vector (i.e. sans singular Γ-irreps)
    nᵀ⁺ᴸs = PBC.sum_symbases.(Ref(sb), cⁱs)
    nᴸ    = sb[idx¹ᴸ]
    global nᵀs   = nᵀ⁺ᴸs .- Ref(nᴸ)
    for nᵀ in nᵀs
        # tricky indexing: this achieves the ordering-adjusted subtraction into the `sb`
        # basis (note the broadcasted assignment, which is what makes this work)
        nᵀ[Γidxs] .-= ms²ᵀ
    end
    unique!(nᵀs) # remove equivalent solutions (due to expansion non-uniqueness)
    io ≠ stdout && println("   ... found ", length(nᵀs), " minimal connectivity solutions")

    # construct human-readable symmetry vectors
    ios = [IOBuffer() for _ in eachindex(nᵀs)]
    prettyprint_symmetryvector.(ios, nᵀs, Ref(sb.irlabs))
    nᵀs_str = String.(take!.(ios))

    # insert (▪)²ᵀ in Γ-irrep spot
    nᵀs_str .= add_content_to_symvec_str_at_kidx.(nᵀs_str, Γidx_in_sb, Ref(repr²ᵀ)) # from /test/text_utils.jl
    if outputtype == "latex" 
        nᵀs_str .= Ref('$'*"\\msf{").*convert_irreplabel2latex.(nᵀs_str).*Ref('}'*'$') # from /test/text_utils.jl
    end

    # write results to screen
    minN   = Nsolutions_found + 1
    plural = length(nᵀs) > 1

    (isℤ₁ && minN == 1) && println(io, "[Trivial symmetry indicator group, ℤ₁]\n")
    println(io,
        "**", minN, ordinal_indicator(minN)," minimal solution:** ",
        "νᵀ = ", νᵀ, " (", length(nᵀs), " symmetry vector", plural ? "s" : "", ")\n")

    flush(io)
    if minN == 1 # only print table for 1st minimal solution
        if !isℤ₁
            # find "Z₂" factor-type topology of each solution
            topos = topology_from_2T1L_xor_1L.(nᵀs, Ref(nᴸ), Ref(ms²ᵀ), Ref(Γidxs), Ref(F))

            io ≠ stdout && println("      ... computed associated xor-topology via EBRs")

            contents = [nᵀs_str topos]
            header   = ["nᵀ", "topology"]
            if all(==(nontrivial), topos)
                println(io, "[Filling-enforced topology]\n")
            end
        else
            contents = nᵀs_str
            header   = ["nᵀ"]
        end
        pretty_table(io, contents, header; # contents & header row
                            alignment = :l,
                            table_config(outputtype)...
                    )
        println(io)
    end
    flush(io)

    # check if we are done
    Nsolutions_found += 1
    Nsolutions_found ≥ Nsolutions && break
    if Nsolutions_found ≥ 1 && (sgnum ∈ (2,10,47))  # skip computational quagmires 
        println(io, "\nSkipped higher-order solutions...")
        break
    end
end # while
println(io)
stdout ≠ io && println("   → total time: ", round((time()-t₀)/60, digits=1), " min")


