"""
    transverse_vrep(lgirs::AbstractVector{LGIrrep{D}}) where D --> LGIrrep{D}

Create an synthetic representation of the virtual representation associated with the
transverse modes at Γ and ω=0. This is a "fake" representation, which will not obey the
representation algebra in general, but which will have the correct irrep dimensionality (2
in 3D, 1 in 2D) and characters. 

The irrep label will be a sum of the constituent irrep labels at Γ and will be distinguished
from other irreps by being enclosed in parentheses.
"""
function transverse_vrep(lgirs::AbstractVector{LGIrrep{D}}) where D
    lg = group(first(lgirs))
    χs = get_symvals²ᵀ(lg) # characters of vrep at ω=0, Γ
    identity_idx = something(findfirst(isone, lg))
    _irdim = convert(Int, χs[identity_idx])
    _irdim == 2 || error("unexpectedly obtained transverse irrep of dimension other than 2")
    diag_entries = χs ./ _irdim
    # create "fake" or "synthetic" representation matrices - which usually will not even
    # realize the algebra of the reps - but which have the correct characters (i.e. trace)
    # at Γ, corresponding to `transverse_vrep_str`; characters are all we need for
    # compatibility analysis
    matrices = [Matrix{Float64}(d*I(_irdim)) for d in diag_entries]
    translations = [zeros(D) for _ in eachindex(lg)]
    cs = something(find_representation(χs, lgirs))
    irlab = '(' * transverse_vrep_str(cs, lgirs) * ')'
    return LGIrrep{D}(irlab, lg, matrices, translations, REAL, false)
end

function transverse_vrep(lgirsd::AbstractDict{String, <:AbstractVector{<:LGIrrep}})
    return transverse_vrep(lgirsd["Γ"])
end

function transverse_vrep_str(cs::AbstractVector{<:Integer}, lgirs::AbstractVector{<:LGIrrep})
    length(cs) == length(lgirs) || error("`cs` and `lgirs` have dissimilar lengths")
    first = true
    io = IOBuffer()
    for (i, c) in enumerate(cs)
        iszero(c) && continue
        sign_str = c < 0 ? "-" : first ? "" : "+"
        first = false
        absc = abs(c)
        c_str = isone(absc) ? "" : string(absc)
        print(io, sign_str, c_str, label(lgirs[i]))
    end
    return String(take!(io))
end

""" 
    is_vrep(lgir::LGIrrep) -> Bool

Return whether or not an irrep is a "fake" virtual representation associated with the
singularity of transverse modes at Γ and ω=0. See [`transverse_vrep`](@ref).
"""
is_vrep(lgir::LGIrrep) = startswith(label(lgir), "(Γ")