using Crystalline
import Crystalline: rotation, rotation_order_3d
using PhotonicBandConnectivity
using PhotonicBandConnectivity: 
        find_representation¹ᴸ, find_representation²ᵀ,
        find_minimum_bandreps_regular2T, 
        find_minimum_bandreps_regular1L
using SymmetryBases
using LinearAlgebra: det, tr
using Test

# ---------------------------------------------------------------------------------------- #

rotation_order(op::SymOperation{3}) = (W=rotation(op); rotation_order_3d(det(W), tr(W)))

# ---------------------------------------------------------------------------------------- #

@testset "Equivalence of regular2T and regular1L approaches" begin 
    for timereversal in (true, false)

        sgnums = collect(1:MAX_SGNUM[3])
        # restrict to sgs that do not have roto-inversions
        filter!(sgnums) do sgnum
            !any(op->rotation_order(op)∈(-1, -3, -4, -6), spacegroup(sgnum, Val(3)))
        end
        # restrict to sgs that have regular 2T representations
        filter!(sgnums) do sgnum
            all(≥(0), find_representation²ᵀ(sgnum, timereversal=timereversal))
        end

        # -------------------------------------------------------------------------------- #

        for sgnum in sgnums
            # irreps at Γ and Hilbert basis
            lgirs = lgirreps(sgnum, Val(3))["Γ"]
            timereversal && (lgirs = realify(lgirs))
            
            # pinned irreps at (Γ,ω=0)
            ms¹ᴸ = find_representation¹ᴸ(lgirs)
            ms²ᵀ = find_representation²ᵀ(lgirs)
            ms   = ms²ᵀ .+ ms¹ᴸ

            # solution via 2T (regular 2T) approach
            cⁱs, ν²ᵀ, sb   = find_minimum_bandreps_regular2T(
                                sgnum, lgirs, timereversal, ms²ᵀ, verbose=false)
            nᵀs  = unique!(sort(sum_symbases.(Ref(sb), cⁱs)))
            
            # solution via 2T+1L & 1L (irregular 2T) approach
            cⁱs′, ν²ᵀ′, _, idx¹ᴸ = find_minimum_bandreps_regular1L(
                                sgnum, lgirs, timereversal, ms¹ᴸ, ms, verbose=false)
            nᵀs′ = unique!(sort(sum_symbases.(Ref(sb), cⁱs′) .- Ref(sb[idx¹ᴸ])))

            # --- testing ---
            @test ν²ᵀ == ν²ᵀ′
            @test nᵀs == nᵀs′
        end

    end # for timereversal
end # testset