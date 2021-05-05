using PhotonicBandConnectivity, Crystalline
using Test

timereversal = true

# test that we obtain the same results when including/excluding non-maximal k-points (e.g.,
# lines, planes)
for sgnum in 1:230
    cⁱs¹, νᵀ¹, _ = minimal_expansion_of_zero_freq_bands(sgnum;
                            timereversal=timereversal, allpaths=true)
    cⁱs², νᵀ², _ = minimal_expansion_of_zero_freq_bands(sgnum;
                            timereversal=timereversal, allpaths=false)

    @test νᵀ¹ == νᵀ²
    @test length(cⁱs¹) == length(cⁱs²)
end
