struct Polynomial{T<:Number}
    p::Vector{T}
end

Base.getindex(p::Polynomial, inds...) = p.p[(inds .- 1)...]

function (::Polynomial{T})(x) where T
    return sum([x^n for n in 0:n] .* p.p)
end

"""
    legendre(n::Integer; length::Integer=0) -> Vector{Float64}

Generates the coefficients of the `n`th Legendre polynomial. The coefficients are ordered such
that the first term is a constant and succeeding terms increase in powers of x.

If the `length` parameter is set to a value greater than n + 1, the output vector will have the
length specified.
"""
function legendre(n::Integer; length::Integer=0)
	# ln = (length == 0) ? n : length
	return sum(
        binomial(n,k) * binomial(n+k, k) * (1/2)^k * [(-1)^(x-k) * binomial(k,x) for x in 0:n]
    for k in 0:n)
end

"""
    legendre(n::Integer, x)

Evaluates the `n`th Legendre polynomial at `x`. Legendre polynomials form part of the description
of the angular nodes in spherical harmonics.
"""
legendre(n::Integer, x) = sum([x^n for n in 0:n] .* legendre(n))

"""
    laguerre(n::Integer; length::Integer=0) -> Vector{Float64}

Generates the coefficients of the `n`th Laguerre polynomial. The coefficients are ordered such
that the first term is a constant and succeeding terms increase in powers of x.

If the `length` parameter is set to a value greater than n + 1, the output vector will have the
length specified.
"""
function laguerre(n::Integer)
    return [(-1)^k * binomial(n,k) * factorial(k) for k in 0:n]
end

"""
    laguerre(n::Integer, x)

Evaluates the `n`th Laguerre polynomial at `x`. Laguerre polynomials form the radial portion of the
wavefunctions of the hydrogen atom.
"""
laguerre(n::Integer, x) = sum([x^n for n in 0:n] .* laguerre(n))

#=
"""
    derivative_matrix(len::integer, n)

Generates the `n`th derivative matrix in a polynomial basis for polynomials of length `len`.
"""
function derivative_matrix(len::Integer, n)
    return [(a == b - 1 ? a : 0) for a in 1:(l+1), b in 1:(l+1)]^n
end
=#