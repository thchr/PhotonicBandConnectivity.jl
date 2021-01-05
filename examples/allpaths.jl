sgnum = 17
timereversal = true
allpaths = true

# solve filling problem
cⁱs, νᵀ, sb, idx¹ᴸ = minimal_expansion_of_zero_freq_bands(sgnum; 
                        timereversal=timereversal, allpaths=allpaths)

lgirsd = get_lgirreps(sgnum, Val(3))
timereversal && (foreach(zip(keys(lgirsd), values(lgirsd))) do (klab, lgirs)
                    lgirsd[klab] = realify(lgirs)
                    end)
lgirs_Γ = lgirsd["Γ"]
Γidxs   = PhotonicBandConnectivity.get_Γidxs(lgirs_Γ, sb)   
ms²ᵀ    = PhotonicBandConnectivity.find_representation²ᵀ(lgirs_Γ)

println(ms²ᵀ)
println(PhotonicBandConnectivity.find_representation²ᵀ⁺¹ᴸ(lgirs_Γ))
#=

nᵀ⁺ᴸs = sum_symbases.(Ref(sb), cⁱs)
nᴸ    = sb[idx¹ᴸ]
nᵀs   = nᵀ⁺ᴸs .- Ref(nᴸ)
for nᵀ in nᵀs
    # get rid of singular irrep content at Γ
    nᵀ[Γidxs] .-= ms²ᵀ
end
irlabs = Crystalline.irreplabels(sb);
νᴸ = sb[idx¹ᴸ][end]

# pretty-print transverse solutions
ios = [IOBuffer() for _ in eachindex(nᵀs)]
Crystalline.prettyprint_symmetryvector.(ios, nᵀs, Ref(irlabs))
nᵀs_str = String.(take!.(ios))

# pretty-print longitudinal solution choice
ioᴸ = IOBuffer()
Crystalline.prettyprint_symmetryvector(ioᴸ, nᴸ, irlabs)
nᴸ_str  = String(take!(ioᴸ))

println("\n", "-"^20, "\nTransverse solutions (νᵀ = $νᵀ)")
println.(Ref("  "), nᵀs_str[1:1])

println("\nLongitudinal solution (νᴸ = $νᴸ)")
println("  ", nᴸ_str)

nothing
=#