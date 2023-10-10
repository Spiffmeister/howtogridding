


"""
    cosinespace(R::Vector{TT},θ::TT,ζ::TT,m::Vector,n::Vector) where TT
Cosine term for cylindrical coordinates
"""
function cosinespace(R::Array{TT},θ::TT,ζ::TT,m::Vector,n::Vector) where TT
    local r = 0.0 :: TT
    for j in eachindex(m)
        for i in eachindex(n)
            r += R[i,j]*cos(convert(TT,m[j])*θ - convert(TT,n[i])*ζ)
        end
    end
    return r
end
function cosinespace(R::Array{TT},θ::AbstractVector{TT},ζ,m,n) where TT
    r = zeros(TT,length(θ))
    for i in eachindex(θ)
        r[i] = cosinespace(R,θ[i],ζ,m,n)
    end
    return r
end

"""
    sinusespace(Z::Vector{TT},θ::TT,ζ::TT,m::Vector,n::Vector) where TT
Sin term for cylindrical coordinates
"""
function sinusespace(Z::Array{TT},θ::TT,ζ::TT,m::Vector{Int},n::Vector{Int}) where TT
    local z = 0.0 :: TT
    for j in eachindex(m)
        for i in eachindex(n)
            z += Z[i,j]*sin(m[j]*θ - n[i]*ζ)
        end
    end
    return z
end
function sinusespace(Z::Array{TT},θ::AbstractVector{TT},ζ,m,n) where TT
    z = zeros(TT,length(θ))
    for i in eachindex(θ)
        z[i] = sinusespace(Z,θ[i],ζ,m,n)
    end
    return z
end




"""
    Torus{TT}
Represents a torus in boundary using Fourier series in `R` and `Z` coordinates as,
``R(θ,ζ) = ∑ᵢⱼ Rᵢⱼ cos(mⱼθ - nᵢζ)``
``Z(θ,ζ) = ∑ᵢⱼ Zᵢⱼ sin(mⱼθ - nᵢζ)``
where `Rᵢⱼ` and `Zᵢⱼ` are the Fourier coefficients and `mⱼ` and `nᵢ` are the Fourier modes.
"""
struct Torus{TT}
    R::Array{TT}
    Z::Array{TT}
    m::Vector{Int}
    n::Vector{Int}
end


function (T::Torus)(θ::TT,ζ::TT) where TT
    return [cosinespace(T.R,θ,ζ,T.m,T.n), sinusespace(T.Z,θ,ζ,T.m,T.n)]
end

