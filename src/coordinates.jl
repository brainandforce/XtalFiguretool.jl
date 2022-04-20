# for Cartesian coordinates, just use SVector{D,T}
# These new types are intended for specific dispatch

"""
    AbstractPolar

Supertype for polar, spherical, or cylindrical coordinates.
"""
abstract type AbstractPolar
end

# Base.getindex(v::AbstractPolar, ind) = getproperty(v, ind)

"""
    r(v::AbstractPolar{T}) -> T

Gets the r component of a polar, spherical, or cylindrical coordinate.
"""
r(v::AbstractPolar) = v.r

"""
    theta(v::AbstractPolar{T}) -> T

Gets the theta component of a polar, spherical, or cylindrical coordinate.
"""
theta(v::AbstractPolar) = v.θ
# TODO: This probably needs to be altered
# LinearAlgebra.norm(v::Union{Polar,Spherical}) = v.r

"""
    Polar{T<:Real}

Represents a polar coordinate in 2D.
"""
struct Polar{T<:Real} <: AbstractPolar
    r::T
    θ::T    # in radians, 0 ≤ θ < 2π
    function Polar(r::Real, θ::Real)
        (r, t) = promote(r, θ % 2π)
        return new{typeof(r)}(r, t)
    end
end

# If only angles are given, assume r is 1
Polar(θ) = Polar(1, θ)
Polar(v::SVector{2,T}) where T<:Real = Polar(norm(v), atan(v[2], v[1]))
convert(::Type{Polar}, v::SVector{2,T}) where T<:Real = Polar(v)

"""
    Spherical{T<:Real}

Represents a spherical coordinate in 3D. The ISO 80000-2:2019 convention is used: θ is the polar
angle (deviation from the z-axis, ranging from 0 to π/2) and ϕ is the azimuthal angle (ranging from
0 to 2π).
"""
struct Spherical{T<:Real} <: AbstractPolar
    r::T
    θ::T    # in radians, 0 ≤ θ < π/2
    ϕ::T    # in radians, 0 ≤ θ < 2π
    function Spherical(r::Real, θ::Real, ϕ::Real)
        (r, t, p) = promote(r, θ % 2π, ϕ % 2π)
        return new{typeof(r)}(r, t, p)
    end
end

Base.abs(v::Spherical) = Spherical(abs(r(v)), theta(v), phi(v))

"""
    phi(v::Spherical{T}) -> T

Gets the phi component of a spherical coordinate.
"""
phi(v::Spherical) = v.ϕ

function Spherical(v::AbstractVector{T}) where T<:Real
    @assert length(v) == 3 "v is not a 3-dimensional vector."
    r = norm(v)
    theta = acos(v[3]/r)
    phi = atan(v[2], v[1])
    return Spherical(r, theta, phi)
end

# If only angles are given, assume r is 1
Spherical(θ, ϕ) = Spherical(1, θ, ϕ)

"""
    cart2sph(v::SVector{3,T}) -> Spherical{T}

Converts a set of 3D Cartesian coordinates into spherical coordinates.
"""
cart2sph(v::AbstractVector) = Spherical(v)
convert(::Type{Spherical}, v::SVector{3,T}) where T<:Real = Spherical(v)

function StaticArrays.SVector(v::Spherical{T}) where T
    x = r(v) * cos(phi(v)) * sin(theta(v))
    y = r(v) * sin(phi(v)) * sin(theta(v))
    z = r(v) * cos(theta(v))
    return SVector{3,T}(x, y, z)
end

"""
    sph2cart(v::Spherical{T}) -> SVector{3,T}

Converts a set of spherical coordinates to Cartesian coordinates.
"""
sph2cart(v::Spherical) = SVector(v)
convert(::Type{SVector}, v::Spherical) = SVector(v)

"""
    sphere(r::Real = 1; grit::Integer = 32) -> Matrix{Spherical{Float64}}

Generates a set of points on a sphere, with radius 1 by default. The points are returned in
spherical coordinates.
"""
function sphere(r::Real = 1; grit::Integer=32)
	theta = range(0, stop=pi, length=grit)
	phi = range(0, stop=2*pi, length=2*grit+1)
	return [Spherical(r, t, p) for t in theta, p in phi]
end
