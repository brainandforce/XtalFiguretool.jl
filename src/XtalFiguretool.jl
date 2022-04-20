module XtalFiguretool

using LinearAlgebra
using StaticArrays
using Xtal

# Primitives for polar and spherical coordinates
include("coordinates.jl")
export AbstractPolar, Polar, Spherical
export r, theta, phi, cart2sph, sph2cart, sphere

# Polynomials often used in chemistry
include("polynomials.jl")

# Needed for plotting spherical harmonics
include("sphericalharmonics.jl")
export Y_real

# Draw methods for 3D objects
include("draw3d.jl")
export Draw3D

include("precompile.jl")

end # module
