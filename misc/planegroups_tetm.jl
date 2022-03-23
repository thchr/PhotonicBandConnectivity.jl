using PhotonicBandConnectivity, SymmetryBases, Crystalline
const PBC = PhotonicBandConnectivity
using LinearAlgebra: det
using PrettyTables

using Crystalline: prettyprint_symmetryvector

# ---------------------------------------------------------------------------------------- #

PLANE2SPACE_NUM = (
           1   #= p1   ⇒ P1   =#, 3   #= p2   ⇒ P2   =#, 6   #= p1m1 ⇒ Pm   =#,
           7   #= p1g1 ⇒ Pc   =#, 8   #= c1m1 ⇒ Cm   =#, 25  #= p2mm ⇒ Pmm2 =#,
           28  #= p2mg ⇒ Pma2 =#, 32  #= p2gg ⇒ Pba2 =#, 35  #= c2mm ⇒ Cmm2 =#,
           75  #= p4   ⇒ P4   =#, 99  #= p4mm ⇒ P4mm =#, 100 #= p4gm ⇒ P4bm =#,
           143 #= p3   ⇒ P3   =#, 156 #= p3m1 ⇒ P3m1 =#, 157 #= p31m ⇒ P31m =#,
           168 #= p6   ⇒ P6   =#, 183 #= p6mm ⇒ P6mm =#,
           )

# ---------------------------------------------------------------------------------------- #

function find_representationᵀᴱᵀᴹ(lgirs, polarization::Symbol)
    lg = group(first(lgirs))
    if polarization == :TE || polarization == :te
        symeigs = det.(rotation.(lg)) # B ∥ ̂z: pseudovector
    elseif polarization == :TM || polarization == :tm
        symeigs = ones(length(lg))    # E ∥ ̂z: vector
    end
    return find_representation(symeigs, lgirs)
end

# ---------------------------------------------------------------------------------------- #

function minimal_expansion_planegroup(sgnum::Integer, timereversal::Bool, polarization::Symbol;
            verbose::Bool=false, allpaths::Bool=false, safetychecks::Bool=false)

    # Irreps at Γ
    lgirs = lgirreps(sgnum, Val(2))["Γ"]
    timereversal && (lgirs = realify(lgirs))

    # Hilbert bases and "default" connectivities
    sb, Γidxs = PBC.compatibility_basis_and_Γidxs(lgirs; timereversal, allpaths)
    νsᴴ       = fillings(sb)
    νᴴₘₐₓ     = maximum(νsᴴ)

    # 2T+1L symmetry constraints
    ms     = find_representationᵀᴱᵀᴹ(lgirs, polarization)
    ntidxs = PBC.find_symmetry_constrained_bases(sb, ms, Γidxs)

    return sb[ntidxs], sb

    throw("Found no valid expansions consistent with constraints")
end

# ---------------------------------------------------------------------------------------- #


##
timereversal = false
table_opts = (tf = tf_unicode, alignment = :l, vlines = :none, hlines = [:begin, 1, :end])

for pgnum in 1:17
    local lgirs
    sgnum = PLANE2SPACE_NUM[pgnum]
    
    lgirs = lgirreps(pgnum, Val(2))["Γ"]
    timereversal && (lgirs = realify(lgirs))
    irlabs = label.(lgirs)
    lg = group(first(lgirs))

    msᵀᴱ = find_representationᵀᴱᵀᴹ(lgirs, :TE)
    msᵀᴹ = find_representationᵀᴱᵀᴹ(lgirs, :TM)
    iridxᵀᴱ = only(findall(≠(0), msᵀᴱ))
    iridxᵀᴹ = only(findall(≠(0), msᵀᴹ))

    nsᵀᴱ, sb = minimal_expansion_planegroup(pgnum, timereversal, :TE)
    nsᵀᴹ, _  = minimal_expansion_planegroup(pgnum, timereversal, :TM)
    
    μsᵀᴱ = last.(nsᵀᴱ)
    μsᵀᴹ = last.(nsᵀᴹ)
    μs = fillings(sb)


    iosᵀᴱ = [IOBuffer() for _ in eachindex(nsᵀᴱ)]
    prettyprint_symmetryvector.(iosᵀᴱ, nsᵀᴱ, Ref(sb.irlabs))
    nsᵀᴱ_str = String.(take!.(iosᵀᴱ))

    iosᵀᴹ = [IOBuffer() for _ in eachindex(nsᵀᴹ)]
    prettyprint_symmetryvector.(iosᵀᴹ, nsᵀᴹ, Ref(sb.irlabs))
    nsᵀᴹ_str = String.(take!.(iosᵀᴹ))

    # calculate topology
    BRS = bandreps(pgnum, 2, timereversal=timereversal)
    B   = matrix(BRS)
    F   = Crystalline.smith(B)

    toposᵀᴱ = calc_detailed_topology.(nsᵀᴱ, Ref(B), Ref(F))
    toposᵀᴹ = calc_detailed_topology.(nsᵀᴹ, Ref(B), Ref(F))


    println("═════ Plane group $pgnum ($(iuc(pgnum, 2)) [parent space group $sgnum] ═════")
    println()

    println("─── TE ───\n", 
            "Pinned Γ-irrep:     ", irlabs[iridxᵀᴱ], "\n",
            "ω=0 connectivities: ", sort!(unique(μsᵀᴱ)))
    pretty_table(stdout,
        [nsᵀᴱ_str μsᵀᴱ toposᵀᴱ],    # contents
        ["nᵀᴱ", "μᵀᴱ", "topology"]; # header row
        table_opts...
    )
    println()

    println("─── TM ───\n", 
            "Pinned Γ-irrep:     ", irlabs[iridxᵀᴹ], "\n",
            "ω=0 connectivities: ", sort!(unique(μsᵀᴹ )))
    pretty_table(stdout,
        [nsᵀᴹ_str μsᵀᴹ toposᵀᴹ],    # contents
        ["nᵀᴹ", "μᵀᴹ", "topology"]; # header row
        table_opts...
    )

    println("\n")

end

