# (sum of) symmetry eigenvalues of ω=0 branches
# in comments below, we refer to a rotation angle θ, which is 2π/n with `n` denoting the 
# rotation order

# two transverse and two longitudinal plane waves (2T+1L)
function get_symval²ᵀ⁺¹ᴸ(
    op::SymOperation{D},
    polarization::Union{Nothing, Symbol} = nothing
) where D
    D == 3 || return get_symval²ᵀ(op, polarization) # no longitudinal modes for D ≤ 2
    W = rotation(op)
    rotval = Crystalline.rotation_order(W)
    n = abs(rotval)
    # This covers every case, including rotations, mirrors, rotoinversions, & inversion
    return sign(rotval) * (2cospi(2/n) + one(Float64))
end
function get_symvals²ᵀ⁺¹ᴸ(
    ops::AbstractVector{<:SymOperation},
    polarization::Union{Nothing, Symbol} = nothing
)
    return get_symval²ᵀ⁺¹ᴸ.(ops, polarization)
end

# single longitudinal plane wave (1L)
get_symval¹ᴸ(::SymOperation{3}) = one(Float64)
get_symval¹ᴸ(::SymOperation) = zero(Float64) # no longitudinal modes for D ≤ 2
get_symvals¹ᴸ(ops::AbstractVector{SymOperation{3}}) = ones(Float64, length(ops))
get_symvals¹ᴸ(ops::AbstractVector{<:SymOperation})  = zeros(Float64, length(ops))

# two transverse plane waves (2T)
function get_symval²ᵀ(op::SymOperation{3}, #=polarization=#::Nothing=nothing)
    W = rotation(op)
    rotval = Crystalline.rotation_order(W)   
    n = abs(rotval) # rotation order 

    if !signbit(rotval)                                     # ← Proper rotation
        # The symmetry eigenvalues are those of of the 2×2 rotation matrix R(θ) ≡ 
        # [c s; c -s] with c ≡ cos(θ), s ≡ sin(θ), and θ ≡ 2π/n, i.e. e⁺ⁱᶿ and e⁻ⁱᶿ
        return 2cospi(2/n) # eⁱᶿ + e⁻ⁱᶿ = 2cos(θ) [θ = 2π/n]

    else                                                    # ← Improper rotation
        # It is not generally possible to infer the all the symmetry eigenvalues of 
        # roto-inversions with rotval = (-1, -3, -4, -6) for the two transverse 
        # plane waves (2T) in isolation. This is because there are there no lines of
        # symmetry from Γ along which 2T could be symmetry-allowed eigenfunctions 
        # for the rotoinversions.
        # Instead, we pick a _possible_ choice for 1L and infer _possible_ symmetry
        # values from [2T+1L] - 1L

        # It _is_ possible for a simple mirror though (i.e., rotation followed by
        # inversion, i.e. -2 === m): the right choice is to pick the symmetry
        # eigenvalues as +1 and -1 (again, assuming two transverse plane waves along
        # each high-symmetry k-vector)
        if rotval == -2                     # ← Mirror
            return zero(Float64) # (+1) + (-1)

        elseif rotval ∈ (-1, -3, -4, -6)    # ← Roto-inversions & inversion
            
            # In general, we have: 
            #   [2T+1L] - 1L = [(-eⁱᶿ) + (-e⁻ⁱᶿ) + (-1)] - (+1) = -2cos(θ) - 2.0
            # For inversion specifically, this is:
            #   [2T+1L] - 1L = [(-1) + (-1) + (-1)] - (+1) = -4
            return -2cospi(2/n) - 2.0
            # SGs w/ non-mirror rotoinversions are 81:88, 111:142, 147:148, 162:167, 
            # 174:176, 187:194, 200:206, and 215:230
        end
    end
end

function get_symval²ᵀ(op::SymOperation{2}, polarization::Symbol)
    if polarization == :TM
        return one(Float64)
    elseif polarization == :TE
        return det(rotation(op))
    else
        throw(DomainError(polarization, "polarization must be :TE or :TM in 2D"))
    end
end

function get_symvals²ᵀ(
    ops::AbstractVector{<:SymOperation},
    polarization::Union{Nothing, Symbol}=nothing
)
    return get_symval²ᵀ.(ops, polarization)
end

# convenience accessors via space/little groups, ensuring primitive basis
for f in (:get_symvals²ᵀ⁺¹ᴸ, :get_symvals²ᵀ)
    @eval function $f(
        sg::Union{LittleGroup, SpaceGroup},
        polarization::Union{Nothing, Symbol} = nothing
    )
        return $f(operations(primitivize(sg)), polarization)
    end
end
get_symvals¹ᴸ(sg::Union{LittleGroup, SpaceGroup}) = get_symvals¹ᴸ(operations(primitivize(sg)))