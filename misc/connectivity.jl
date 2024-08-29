using Crystalline
using Graphs, MetaGraphs
using LinearAlgebra: pinv
using StaticArrays

#=
    Many of the methods in this file are generalizations of similar very similar methods in
    Crystalline/src/compatibility.jl.
    The purpose of this file is to allow methods to give an indication of the connectedness
    between high-symmetry points in the Brillouin zone, and how they are connected. 
    j
    TODO: These methods should probably eventually be merged with the functionality in 
          Crystalline/src/compatibility.jl. At present, it is just easier to have it 
          separated, because the use-case here is quite specific to analyzing the validity
          of the symmetry vectors extracted by PhotonicBandConnectivity
=#

# ---------------------------------------------------------------------------------------- #
# METHODS

function is_compatible(kv::KVec{D}, kv′::KVec{D}, cntr::Char) where D
    # This method determines whether there is a solution to 
    #   kv + G = kv′(αβγ′)
    # and if so, returns G and αβγ′; otherwise, returns `nothing`. `kv` must be special.
    isspecial(kv′) && return nothing # must have a free parameter to match a special point `kv`
    isspecial(kv) || throw(DomainError(kv, "must be special"))

    k₀, _      = parts(primitivize(kv,  cntr)) 
    k₀′, kabc′ = parts(primitivize(kv′, cntr))

    # check if there's a solution of `kabc′*αβγ′ = k₀ - k₀′` modulo `G` in `αβγ′`
    for Gtup in Iterators.product(ntuple(Returns((0,1,-1)), Val(D))...) # modulo G checks...
        G = SVector{D,Int}(Gtup)
        k₀G = k₀ + G
        αβγ′ = pinv(kabc′)*(k₀G-k₀′) # least squares solution
        k′ = k₀′ + kabc′*αβγ′
        # check if least squares solution actually is a solution
        compat_bool = isapprox(k₀G, k′, atol=Crystalline.DEFAULT_ATOL)
        if compat_bool
            return αβγ′, G
        end
    end

    return nothing
end

function find_compatible(kv::KVec{D}, kvs′::AbstractVector{KVec{D}}, cntr::Char) where D
    !isspecial(kv) && throw(DomainError(kv, "input kv must be a special k-point"))

    compat_idxs = Int64[]
    compat_αβγs = SVector{D,Float64}[]
    compat_Gs   = SVector{D,Int}[]
    @inbounds for (idx′, kv′) in enumerate(kvs′)
        αβγ′_and_G = is_compatible(kv, kv′, cntr)
        if !isnothing(αβγ′_and_G)
            αβγ′, G = αβγ′_and_G
            push!(compat_idxs, idx′)
            push!(compat_αβγs, αβγ′)
            push!(compat_Gs, G)
        end
    end

    return compat_idxs, compat_αβγs, compat_Gs
end

# let kvsᴬ refer to a set of special k-points and kvsᴮ to a set of mixed k-points (special
# and non-special)
function connectivity((kvsᴬ, klabsᴬ), (kvsᴮ, klabsᴮ), cntr)

    all(isspecial, kvsᴬ) || throw(DomainError(kvsᴬ, "must only include special points"))

    Nkᴬ = length(kvsᴬ)
    D = dim(first(kvsᴬ))

    # find compatible vectors between A and B
    compat_idxs = Vector{Vector{Int}}(undef, Nkᴬ)
    compat_αβγs = Vector{Vector{SVector{D,Float64}}}(undef, Nkᴬ)
    compat_Gs = Vector{Vector{SVector{D,Int}}}(undef, Nkᴬ)
    cgraph = MetaGraph(Nkᴬ) # connectivity graph for special k-vecs
    for (idxᴬ, (kvᴬ, klabᴬ)) in enumerate(zip(kvsᴬ, klabsᴬ))
        compat_idxs[idxᴬ], compat_αβγs[idxᴬ], compat_Gs[idxᴬ] = find_compatible(kvᴬ, kvsᴮ, cntr)
        set_props!(cgraph, idxᴬ, Dict(:kv => kvᴬ, :klab => klabᴬ)) # TODO: add idx of kvᴬ in kvsᴮ?
    end

    for (idxᴬ¹, kvᴬ¹) in enumerate(kvsᴬ)
        compat_idxs¹ = compat_idxs[idxᴬ¹]
        for (idxᴬ², kvᴬ²) in enumerate(kvsᴬ)
            idxᴬ² ≤ idxᴬ¹ && continue # avoid double & self-comparisons
            compat_idxs² = compat_idxs[idxᴬ²]

            # find overlap of compat_idxs¹ and compat_idxs² which are lines
            local_idxsᴮ = findall(compat_idxs²) do idxᴮ
                (idxᴮ∈compat_idxs¹) & (Crystalline.nfreeparams(kvsᴮ[idxᴮ]) == 1)
            end
            isempty(local_idxsᴮ) && continue
            idxsᴮ = getindex.(Ref(compat_idxs²), local_idxsᴮ)

            for idxᴮ in idxsᴮ
                println(
                    klabsᴬ[idxᴬ¹], " = ", string(kvsᴬ[idxᴬ¹]), " connected to ",
                    klabsᴬ[idxᴬ²], " = ", string(kvsᴬ[idxᴬ²]), " via ",
                    klabsᴮ[idxᴮ], " = ", string(kvsᴮ[idxᴮ]))
            end
            add_edge!(cgraph, idxᴬ¹, idxᴬ²)
            set_props!(cgraph, Edge(idxᴬ¹, idxᴬ²), 
                Dict(:klabs => getindex.(Ref(klabsᴮ), idxsᴮ),
                     :kvs   => getindex.(Ref(kvsᴮ), idxsᴮ),
                     :kidxs => idxsᴮ)
               )
        end
    end
   
    return cgraph
end

## --------------------------------------------------------------------------------------- #
# SCRIPTING
timereversal = true
sgnum = 230

brs = bandreps(sgnum, 3; timereversal)
lgirsd = lgirreps(sgnum, Val(3))
if timereversal
    lgirsd = Dict(klab => realify(lgirs) for (klab, lgirs) in lgirsd)
end
kvsᴬ, klabsᴬ = brs.kvs, brs.klabs
kvsᴮ, klabsᴮ = position.(first.(values(lgirsd))), klabel.(first.(values(lgirsd)))
cntr = centering(sgnum)
cg = connectivity((kvsᴬ, klabsᴬ), (kvsᴮ, klabsᴮ), cntr)

# ---------------------------------------------------------------------------------------- #
# PLOT GRAPH
using Plots, GraphRecipes, NetworkLayout


GraphRecipes.graphplot(cg;
    names = Dict(k=>v[:klab] for (k,v) in cg.vprops),
    edgelabel = Dict(Pair(e)=>join(v[:klabs],", ") for (e,v) in cg.eprops),
    edgelabel_offset=.0,
    nodeshape=:circle,
    nodesize=0.15,
    curvature_scale=.1,#curves=false,
    method=:spring,
    method_kw=Dict(:initialpos => [ones(2) for _ in cg.vprops])
    )

# ---------------------------------------------------------------------------------------- #
# PLOT GRAPH
#using GraphMakie, GLMakie
#GraphMakie.graphplot(cg;
#    nlabels = [cg.vprops[v][:klab] for v in 1:length(cg.vprops)],
#    nlabels_align=(:center,:center),
#    node_size=35,
#    node_color=GLMakie.Colors.colorant"gray",
#    elabels = [join(cg.eprops[e][:klabs], ", ") for e in edges(cg)]
#    )