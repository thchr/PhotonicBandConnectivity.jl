using Crystalline, PhotonicBandConnectivity

# compute all symmetry vectors of photonic bands intrinsically connected to ω=0
D = 3
timereversal = true
nsv = Vector{Vector{SymmetryVector{D}}}(undef, 230)
for sgnum in 1:230
    t₀ = time()
    print("#", sgnum)
    ns = transverse_symmetry_vectors(sgnum, Val(D); timereversal)
    nsv[sgnum] = ns
    println(" (", round(time()-t₀, digits=1), " s): ", length(ns), " solutions (",
            minimum(occupation, ns), " ≤ μ ≤ ", maximum(occupation, ns), ")")
end
nsv