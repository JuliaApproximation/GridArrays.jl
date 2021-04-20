
"""
A `ProductGrid` represents the cartesian product of other grids.

`struct ProductGrid{TG,T,N} <: AbstractGrid{T,N}`

Parameters:
- TG is a tuple of (grid) types
- T is the element type of the grid
- N is the dimension of the grid layout
"""
struct ProductGrid{TG,T,N} <: AbstractGrid{T,N}
	grids	::	TG
end

# Generic functions for composite types:
components(grid::ProductGrid) = grid.grids
# component(grid::ProductGrid, range::AbstractRange) = productgrid(grid.grids[range]...)
component(grid::ProductGrid, range::AbstractRange) = error("Deprecated, make a product grid instead.")

function ProductGrid(grids...)
	T = mapreduce(numtype, promote_type, grids)
	N = sum(map(dimension, grids))
	ProductGrid{typeof(grids),SVector{N,T},length(grids)}(grids)
end

size(g::ProductGrid) = map(length, g.grids)
size(g::ProductGrid, j::Int) = length(g.grids[j])

covering(g::ProductGrid) = productdomain(map(covering, components(g))...)
isperiodic(g::ProductGrid) = reduce(&, map(isperiodic, components(g)))

function unsafe_grid_getindex(g::ProductGrid{TG,T,N}, I::Vararg{Int,N}) where {TG,T,N}
	@inbounds convert(T, map(getindex, g.grids, I))
end

similargrid(grid::ProductGrid, ::Type{T}, dims...) where T = error()#ProductGrid([similargrid(g, eltype(T), dims[i]) for (i,g) in enumerate(components(grid))]...)

toiterator(::Type{ProductGrid}, a::AbstractGrid) = (a,)
toiterator(::Type{ProductGrid}, a::ProductGrid) = Iterators.flatten(ntuple(k->toiterator(ProductGrid,component(a,k)),Val(dimension(a))))
toiterator(::Type{ProductGrid}, components::Vararg{AbstractGrid,N}) where N = Iterators.flatten(ntuple(k->toiterator(ProductGrid,components[k]),Val(N)))
flatten(::Type{ProductGrid}, components::Vararg{AbstractGrid,N}) where N = flatten(ProductGrid, toiterator(ProductGrid, components...))
function flatten(::Type{ProductGrid}, a::Iterators.Flatten)
	collect(a)
end

cross(grids::AbstractGrid...) = productgrid(grids...)
^(grid::AbstractGrid, n::Int) = productgrid(grid, n)
# In order to avoid strange nested structures, we flatten the arguments
productgrid(grids::AbstractGrid...) = ProductGrid(flatten(ProductGrid, grids...)...)
productgrid(grid::AbstractGrid, n::Int) = ProductGrid(ntuple(k->grid,Val(n))...)
# Disallow cartesian products with just one argument
productgrid(grid::AbstractGrid) = grid
