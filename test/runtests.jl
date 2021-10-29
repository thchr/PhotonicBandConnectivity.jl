using Crystalline
using PhotonicBandConnectivity
using Test

@testset "Pinned irreps for 2T, 1L, and 2T+1L" begin 
    for sgnum in 1:230
        for timereversal in (true, false)
            lgirs = lgirreps(sgnum, Val(3))["Γ"]
            timereversal && (lgirs = realify(lgirs))

            # pinned irreps at (Γ,ω=0)
            ms¹ᴸ = PhotonicBandConnectivity.find_representation¹ᴸ(lgirs)
            ms²ᵀ = PhotonicBandConnectivity.find_representation²ᵀ(lgirs)
            ms   = PhotonicBandConnectivity.find_representation²ᵀ⁺¹ᴸ(lgirs)
            
            @test all(≥(0), ms)             # check: 2T+1L irreps regular
            @test all(≥(0), ms¹ᴸ)           # check: 1L irrep regular (Γ₁)
            @test all(ms .== ms²ᵀ .+ ms¹ᴸ)  # check: [2T+1L] = [2T] + [1L]
        end
    end
end

# other tests
include("check_pick1L_invariance.jl")
include("regular2T_vs_regular1L_equivalence.jl")