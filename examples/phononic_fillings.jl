using Pkg
dirname(Pkg.project().path) == (@__DIR__) || (Pkg.activate(@__DIR__); cd(@__DIR__))

using Crystalline
using Crystalline: prettyprint_symmetryvector
using SymmetryBases
using PhotonicBandConnectivity
using PhotonicBandConnectivity: minimal_expansion_of_zero_freq²ᵀ⁺¹ᴸ_bands
using Nemo: MatrixSpace, ZZ
using ProgressMeter
using JLD2

# find all the minimum fillings below the fundamental gap in systems where the longitudinal
# bands are real, such as in phononic or magnonic systems
sgnums   = 1:230
has_tr   = true
workhard = true
savedata = false

μ²ᵀ⁺¹ᴸs = Vector{Int}(undef, length(sgnums)) # minimal fillings below fundamental gap
μs      = Vector{Int}(undef, length(sgnums))
topos   = Vector{Vector{TopologyKind}}(undef, length(sgnums))
n²ᵀ⁺¹ᴸs_str = Vector{Vector{String}}(undef, length(sgnums))
for (sgidx, sgnum) in enumerate(sgnums)
    print("SG ", sgnum, ":\t")
    
    cⁱs, μ²ᵀ⁺¹ᴸs[sgidx], sb = minimal_expansion_of_zero_freq²ᵀ⁺¹ᴸ_bands(sgnum, has_tr)
    μs[sgidx] = minimum(fillings(sb)) # "default" fillings
    print("μ²ᵀ⁺¹ᴸ = ", μ²ᵀ⁺¹ᴸs[sgidx], " (μ = ", μs[sgidx], ")")

    # calculate unique symmetry vectors
    n²ᵀ⁺¹ᴸs = unique!(sort(sum_symbases.(Ref(sb), cⁱs)))
    print(" w/ ", length(n²ᵀ⁺¹ᴸs), " solution", length(n²ᵀ⁺¹ᴸs)≠1 ? "s" : "")

    # pretty-print the symmetry vector solutions
    ios = [IOBuffer() for _ in eachindex(n²ᵀ⁺¹ᴸs)]
    prettyprint_symmetryvector.(ios, n²ᵀ⁺¹ᴸs, Ref(sb.irlabs))
    n²ᵀ⁺¹ᴸs_str[sgidx] = String.(take!.(ios))

    # bail out on computational quagmires for topology evaluation
    if (!workhard &&
        # the tr-broken specific ones are just extremely slow to compute nontopo_sb for; 
        # didn't check fragile phases for those (didn't terminate in +24h for any case)
        (has_tr && sgnum ∈ (2,10,47) || (!has_tr && sgnum ∈ (83,174,175,176))) )

        topos[sgidx] = TopologyKind[]
        println("\t [skipped...]"); continue
    end
    
    # calculate topology
    BRS = bandreps(sgnum, 3, timereversal=has_tr)
    B   = matrix(BRS, true)
    F   = Crystalline.smith(B)
    topo_class = classification(BRS)
    print("\t [", topo_class)

    nontopo_sb = nontopological_basis(F, BRS)
    trivial_idxs, fragile_idxs = split_fragiletrivial(nontopo_sb, B)
    println(!isempty(fragile_idxs) ? "+fragile" : "", "]")

    flush(stdout)
    if isempty(fragile_idxs)
        topo_class == "Z₁" && (topos[sgidx] = fill(trivial, length(n²ᵀ⁺¹ᴸs)); continue)

        Bℤ       = MatrixSpace(ZZ, size(B)...)(B)  # Nemo-version of matrix
        topos[sgidx] = @showprogress map(n -> calc_topology(n, Bℤ), n²ᵀ⁺¹ᴸs)
    else
        M         = matrix(sb)
        nontopo_M = matrix(nontopo_sb)
        trivial_M = nontopo_M[:,trivial_idxs]
        topos[sgidx] = @showprogress map(n²ᵀ⁺¹ᴸs) do n
                            calc_detailed_topology(n, nontopo_M, trivial_M, 
                                                   workhard ? M : nothing)
                       end
    end

    # print if filling-enforced topology detected
    if all(==(nontrivial), topos[sgidx])
        println("\t⇒ Filling-enforced stable topology!")
    elseif all(==(fragile), topos[sgidx])
        println("\t⇒ Filling-enforced fragile topology!")
    elseif all(≠(trivial), topos[sgidx])
        println("\t⇒ Filling-enforced mixed stable & fragile topology")
    end
end

# save data, if requested
if savedata && sgnums == 1:230 && has_tr
    @save "phonon-connectivity-result-with-tr.jld2" sgnums has_tr workhard μ²ᵀ⁺¹ᴸs μs topos
elseif savedata && sgnums == 1:230 && !has_tr
    @save "phonon-connectivity-result-without-tr.jld2" sgnums has_tr workhard μ²ᵀ⁺¹ᴸs μs topos
end