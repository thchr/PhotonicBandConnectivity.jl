using Pkg
dirname(Pkg.project().path) == (@__DIR__) || (Pkg.activate(@__DIR__); cd(@__DIR__))

using Crystalline
using SymmetryBases
using PhotonicBandConnectivity
using PhotonicBandConnectivity: minimal_expansion_of_zero_freq²ᵀ⁺¹ᴸ_bands
using Nemo: MatrixSpace, ZZ
using ProgressMeter
using JLD2

# pretty much the same purpose as `examples/phononic_fillings.jl`, except we here focus on
# the computational quagmires in the tr-broken case for space groups 83, 174, 175, and 176. 
# for these space groups, computing `nontopological_bases(...)` is pretty much impossible
# (at least, doesn't terminate in +24 h on a laptop), so to verify that they don't host 
# filling enforced topology, we just check that they indeed have both trivial and nontrivial
# states in their solution space. We do this by using `calc_detailed_topology` but with 
# input `nontopo_M` simply substituted for the EBR basis `M`, instead of the actual Hilbert
# basis for all nontopological states. Doing it this way means that the output is somewhat
# misleading: `trivial` classification still refers to trivial, but `nontrivial` could refer
# to either genuinely (i.e. stable) _or_ fragile nontriviality (i.e. `nontrivial` now may
# effectively also mean `fragile`).
# still, this is fine for us: we just need to verify that the solution space contains at
# least one trivial state, so that we can rule out that they contain filling-enforced 
# topology: indeed, this is what we find when we run this script.


sgnums   = (83, 174, 175, 176)
has_tr   = false

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
 
    # get EBRS
    BRS = bandreps(sgnum, 3, timereversal=has_tr)
    B   = matrix(BRS, true)
    F   = Crystalline.smith(B)
    topo_class = classification(BRS)
    println("\t [", topo_class, "(+maybe fragility; unchecked)]")
    flush(stdout)

    # check for presence of trivial states (bit "faux", but kosher for our needs; see
    # introductory text in this file)
    nontopo_M = B  # pretend there's no possibility for fragile phases
    trivial_M = nothing
    M         = matrix(sb)
    topos[sgidx] = @showprogress map(n²ᵀ⁺¹ᴸs) do n
                            calc_detailed_topology(n, nontopo_M, trivial_M, M)
                   end

    # print if filling-enforced topology detected
    if all(==(nontrivial), topos[sgidx])
        println("\t⇒ Filling-enforced topology (may include fragile phases)!")
    end
end