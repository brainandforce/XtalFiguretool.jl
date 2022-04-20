
"""
    Pcos(l::Integer, m::Integer, theta)

Evaluates the associated Legendre polynomial of `cos(theta)`, corresponding to one of the terms
needed to describe an orbital with angular momentum `l` and magnetic quantum number `m`.
"""
function Pcos(l::Integer, m::Integer, theta)
    # Derivative matrix
    deriv = [(a == b - 1 ? a : 0) for a in 1:(l+1), b in 1:(l+1)]
    # Condon–Shortley phase is applied HERE
    poly = (-1)^m * (deriv^m * legendre(l))[1:(end-m)]
    return sum(poly .* [cos(theta)^n for n = 0:(l-m)]) * sin(theta)^m
end

"""
    Y_real(l::Integer, m::Integer, v)

Evaluates the real spherical harmonic at the 3-dimensional vector `v`.
"""
function Y_real(l::Integer, m::Integer, v::Spherical)
    mabs = abs(m)
    @assert mabs <= l "|m| must be greater than or equal to l"
    rho = r(v) * sqrt((2l+1)/(4pi)) * Pcos(l, mabs, theta(v))
    # This can return immediately if m == 0
    m == 0 && return rho
    # Condon–Shortley phase is applied in `Pcos()`
    c = sqrt(2*factorial(l - mabs)/factorial(l + mabs))
    return rho * c * (m > 0 ? cos(mabs * phi(v)) : sin(mabs * phi(v)))
end

Y_real(l::Integer, m::Integer, v::AbstractVector{T}) where T<:Real = Y_real(l, m, cart2sph(v))

"""
    Y_real(l::Integer, m::Integer, r::Real=1; grit=32) -> Matrix{Spherical{Float64}}

Generates a set of points that can be used to plot a spherical harmonic.
"""
function Y_real(l::Integer, m::Integer, r::Real = 1; grit::Integer=32)
    # Generate the sphere that has all the data points
    a = sphere(r; grit)
    # Calculate the spherical harmonic radial values throughout the sphere
    rho = Y_real.(l, m, a)
    # Convert to spherical coordinates
    return Spherical.(rho, theta.(a), phi.(a))
end

"""
    Y_real(s::SphericalComponents; grit::Integer = 32) = Matrix{Spherical{Float64}}

Evaluates a linear combination of spherical harmonics on the surface of a sphere generated with
`sphere()`.
"""
function Y_real(s::SphericalComponents{Lmax}; grit::Integer = 32) where Lmax
    sph = sphere(;grit)
    # Preallocate the output matrix
    r = zeros(size(sph))
    for l in 0:Lmax
        for m in -l:l
            r += s[l,m] * Y_real.(l, m, sph)
        end
    end
    return Spherical.(r, theta.(sph), phi.(sph))
end
