
"""
    abstract type SubGrid{T,N} <: AbstractGrid{T,N} end

A subgrid of an underlying grid.
"""
abstract type SubGrid{T,N} <: SimpleLazyGrid{T,N} end

supergrid(g::SubGrid) = g.supergrid

include("indexed.jl")
include("masked.jl")
include("product.jl")
include("boundary.jl")


subgrid(grid::AbstractGrid, domain::Domain) = MaskedGrid(grid, domain)

function subgrid(grid::AbstractEquispacedGrid, domain::AbstractInterval)
    a = infimum(domain)
    b = supremum(domain)
    h = step(grid)
    idx_a = convert(Int, ceil( (a-grid[1])/step(grid))+1 )
    idx_b = convert(Int, floor( (b-grid[1])/step(grid))+1 )
    idx_a = max(idx_a, 1)
    idx_b = min(idx_b, length(grid))
    IndexSubGrid(grid, idx_a:idx_b, domain)
end

function subgrid(grid::ScatteredGrid, domain::Domain)
    mask = in.(grid, Ref(domain))
    points = grid.points[mask]
    ScatteredGrid(points, domain)
end

function subgrid(grid::ProductGrid, domain::ProductDomain)
    if ncomponents(grid) == ncomponents(domain)
        productgrid(map(subgrid, components(grid), components(domain))...)
    else
        MaskedGrid(grid, domain)
    end
end
