using PhotonicBandConnectivity
using Test
include((@__DIR__)*"/../misc/watanabelu_results.jl") # loads Watanabe & Lu data (in `Msᵂᴸ`)
using Main.WatanabeLuResults

sgnums = 1:230
data   = minimal_expansion_of_zero_freq_bands.(sgnums, timereversal=true, verbose=true)
νᵀs    = getindex.(data, 2)

# Compare with Watanabe & Lu
Q = [[sg, M, Mbound] for (sg, M, Mbound) ∈ zip(sgnums, νᵀs, Msᵂᴸ)]
issues      = map(x-> x[2]≥(x[3]) ? " " : "!",  Q)
differences = map(x-> x[2]==(x[3]) ? " " : "≠", Q)

foreach(vcat.(Q, issues, differences)) do x
    println(
        "|", " "^(4-ndigits(x[1])), x[1], " |", " "^(3-ndigits(x[2])),  # sgnum & spacing
        x[2], " | ",                        # our M predictions
        x[3] == 2 ? "=" : "≥", x[3], " | ", # M-bounds from Watanabe & Lu
        x[4], " | ",                        # bound violations
        x[5], " |"                          # differences from W&L bound?
    )
end