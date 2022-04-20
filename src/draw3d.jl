"""
    Draw3D

A submodule containing methods for plotting 3D objects.
"""
module Draw3D

using ..XtalFiguretool
using StaticArrays
using Colors
using Xtal
import WGLMakie

"""
    Draw3D.spherical_harmonic(
        s::SphericalComponents;
        offset::AbstractVector{<:Real} = [0, 0, 0],
        grit = 32,
        gradient = false,
        kwargs...
    ) -> Makie.FigureAxisPlot

Creates a surface corresponding to a linear combination of real spherical harmonics.
"""
function spherical_harmonic(
    s::SphericalComponents{Lmax};
    offset::AbstractVector{<:Real} = [0, 0, 0],
    grit = 32,
    gradient = true,
    kwargs...
) where Lmax
    sph = Y_real(s; grit)
    cart = sph2cart.(abs.(sph))
    (x,y,z) = ([v[n] for v in cart] .+ offset[n] for n in 1:3)
    # Automatically determine maximum and minimum for gradient
    if gradient
        peak = maximum(abs.(r.(sph)))
        crange = (-peak, peak)
    else
        crange = (-1, 1)
    end
    return WGLMakie.surface(
        x, y, z,
        color = (gradient ? r.(sph) : sign.(r.(sph))),
        colormap = :bwr,
        colorrange = crange,
        shading = false,
        figure = (; resolution = (1000, 900));
        kwargs...
    )
end

"""
    Draw3D.cp_lobe(s::SphericalComponents; kwargs) -> Makie.FigureAxisPlot

Draws the projection of chemical pressures from `CPpackage2` onto spherical harmonics. This is
identical to `Draw3D.spherical_harmonic()` but uses visualization defaults that are common to
chemical pressure schemes (e.g. white positive and black negative lobes).
"""
cp_lobe(s::SphericalComponents; kwargs...) = spherical_harmonic(s, colormap = :grays; kwargs...)

"""
    Draw3D.unit_cell(b::BasisVectors{3}; kwargs...)

Draws the unit cell of a crystal with basis vectors `b`.
"""
function unit_cell(
    bv::BasisVectors{3};
    color_axes = false,
    kwargs...
)
    # Generate the positions
    reduced_coords = [m .% 2 .^ (1:3) .>= 2 .^ (0:2) for m in 0:7]
    coords = bv .* reduced_coords
    (x,y,z) = ([v[n] for v in coords] for n in 1:3)
    return WGLMakie.wireframe(
        x, y, z;
        kwargs...
    )
end

"""
    Draw3D.scalar_datagrid(g::RealSpaceDataGrid{3,T<:Real}; kwargs...) -> Makie.FigureAxisPlot

Draws a volume plot for a scalar datagrid.
"""
function scalar_datagrid(
    g::RealSpaceDataGrid{3,T};
    kwargs...
) where T<:Real
    b = basis(g)
    gsz = gridsize(g)
    # Get the x,y,z positions in the unit cell
    positions = [
        (ind .- 1) ./ gsz * (basis[n] for n in 1:3)
        for ind in Iterators.product([1:d for d in gsz])
    ]
    (x,y,z) = ([v[n] for v in positions] for n in 1:3)
    return WGLMakie.volume(
        x, y, z,
        shading = false;
        kwargs...
    )
end

end
