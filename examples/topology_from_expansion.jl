using Crystalline, LinearAlgebra, Test
using Crystalline: rotation, rotation_order_3d,
                   prettyprint_symmetryvector
using PhotonicBandConnectivity
using SymmetryBases
using PrettyTables: pretty_table
include("text_utils.jl"); using Main.TextUtils

#=  
    This file analyses the trivial/fragile/nontrivial topology of those spacegroups (68 in 
    total) that have well-defined and regular Γ-irreps at ω=0, such that their symmetry 
    vectors are regular as well. In that case, we can assess the topology of the symmetry
    vectors simply by checking against the Hilbert bases, using `calc_detailed_topology`
    from SymmetryBases.
    For a more general approach, see examples/topology_from_xor.jl, which doesn't require
    regular symmetry content at ω=0 (but cannot treat fragile topology).
=#

has_tr = true
rotation_order(op::SymOperation{3}) = (W=rotation(op); rotation_order_3d(det(W), tr(W)))

# --- find SGs with "regular" ω=0-connected Γ-irrep content ---
sgnums = collect(1:MAX_SGNUM[3])
# restrict to sgs that do not have roto-inversions
filter!(sgnums) do sgnum
    !any(op->rotation_order(op)∈(-1, -3, -4, -6), spacegroup(sgnum, Val(3)))
end
# restrict to sgs that have regular 2T representations
filter!(sgnums) do sgnum
    all(≥(0), PhotonicBandConnectivity.find_representation²ᵀ(sgnum, timereversal=has_tr))
end

# for each of the band-solution in `nᵀs`, check whether or not it can be expanded solely 
# using the nontopological basis (then it is a trivial band-combination!); if not, it must
# be a nontrivial band-combination; also check for fragile topology
filepath = "/mnt/c/Dropbox (MIT)/Web/thchr.github.io/_assets/fragile_connectivity_tables.md" # specify output path
io = open(filepath, "w+")
for (sgidx, sgnum) in enumerate(sgnums)
    # --- analyze topology of minimal connectivities ---

    # band representations
    brs = bandreps(sgnum, 3; timereversal=has_tr)
    B   = stack(brs)
    F   = smith(B)

    # get Hilbert bases and minimal connectivities + symmetry vectors
    cⁱs, νᵀ, sb, _ = minimal_expansion_of_zero_freq_bands.(sgnum, timereversal=has_tr, verbose=false)
    nontopo_sb, _  = nontopological_basis(F, brs)

    # compute the unique symmetry vectors for each compatibility-constrained band solution
    nᵀs = unique!(sort(sum_symbases.(Ref(sb), cⁱs)))

    # determine the trivial, respectively fragile indices in nontopo_sb
    trivial_idxs, fragile_idxs = split_fragiletrivial(nontopo_sb, brs)
    can_be_fragile = !isempty(fragile_idxs)

    # stop early if this SG cannot host any phases fragile at all anyway
    can_be_fragile || continue

    # compute topology of each solution
    topos  = calc_detailed_topology.(nᵀs, Ref(B), Ref(F))

    # --- print results as table ---
    println(io, "## SG ", sgnum, "\n\nνᵀ = ", νᵀ, "\n")

    # print global properties of the space group wrt. topology
        #classification(brs) ≠ "Z₁" && print(" [+NONTRIVIAL]")
        #!isempty(fragile_idxs)     && print(" [+FRAGILE]")

    # construct human-readable symmetry vectors, with the ω=0-connected Γ-irreps highlighted

    # irreps at Γ
    lgirs_Γ = lgirreps(sgnum, Val(3))["Γ"]
    has_tr && (lgirs_Γ = realify(lgirs_Γ))
    irlabs_Γ = label.(lgirs_Γ)
    Γidxs = PhotonicBandConnectivity.get_Γidxs(lgirs_Γ, sb) # mapping between sb and lgirs_Γ

    # string-representation of ω=0 Γ-irreps (in `ms²ᵀ`)
    ms²ᵀ    = PhotonicBandConnectivity.find_representation²ᵀ(lgirs_Γ) # (sorted w/ irlabs_Γ)
    io_ms²ᵀ = IOBuffer()
    prettyprint_symmetryvector(io_ms²ᵀ, ms²ᵀ, irlabs_Γ)
    ms²ᵀ_str = strip(String(take!(io_ms²ᵀ)), ('[', ']'))

    # build a symmetry vector with highlighted ω=0 Γ-irreps: 
    # ... first, subtract ω=0 Γ-irreps from symmetry vector
    nᵀs′    = deepcopy(nᵀs)    
    for nᵀ′ in nᵀs′ 
        # tricky indexing: this achieves the ordering-adjusted subtraction into the `sb`
        # basis (note the broadcasted assignment, which is what makes this work)
        nᵀ′[Γidxs] .-= ms²ᵀ
    end
    # ... now, pretty-print nᵀs′ and add ω=0 Γ-irreps back in with a highlight
    Γidx_in_sb = findfirst(==("Γ"), sb.klabs)
    ios = [IOBuffer() for _ in eachindex(nᵀs′)]
    prettyprint_symmetryvector.(ios, nᵀs′, Ref(irreplabels(sb)))
    nᵀs_str = String.(take!.(ios))
    nᵀs_str .= add_content_to_symvec_str_at_kidx.(nᵀs_str, Γidx_in_sb, Ref("("*ms²ᵀ_str*")²ᵀ"))

    # print as table
    pretty_table(io,
        [nᵀs_str topos],    # contents
        ["nᵀ", "topology"]; # header row
        #tf = tf_unicode,
        #vlines = :none, hlines = [:begin, 1, :end],
        alignment = :l,
        crop = :none,
        tf = tf_markdown,
    )
    println(io)
end
io ≠ stdout && close(io)

# some manual clean-up
run(`sed -i 's/|-/|:/g' $filepath`)
run(`sed -i 's/²ᵀ/\\T/g' $filepath`)

nothing