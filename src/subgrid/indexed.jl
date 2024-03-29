"""
	struct IndexSubGrid{G,I,T,N} <: SubGrid{T,N}

An IndexSubGrid is a subgrid corresponding to a certain range of indices of the
underlying grid.
It is assumed to be an 1D grid.
"""
struct IndexSubGrid{G,I,T,N} <: SubGrid{T,N}
	supergrid  :: G
	subindices :: I
	domain 	   :: Domain

	function IndexSubGrid{G,I,T,N}(supergrid::AbstractArray{T,N}, subindices, domain=Interval(first(supergrid), last(supergrid))) where {G,I,T,N}
		@assert length(subindices) <= length(supergrid)
		new(supergrid, subindices, domain)
	end
end

IndexSubGrid(grid::AbstractArray{T,N}, I, domain=Interval(first(grid), last(grid))) where {T,N} =
    IndexSubGrid{typeof(grid),typeof(I),T,N}(grid, I, domain)

subindices(g::IndexSubGrid) = g.subindices

similar_subgrid(g::IndexSubGrid, g2::AbstractArray) = IndexSubGrid(g2, subindices(g))

length(g::IndexSubGrid) = length(subindices(g))

size(g::IndexSubGrid) = (length(g),)

eachindex(g::IndexSubGrid) = eachindex(subindices(g))

# The speed of this routine is the main reason why supergrid and subindices
# are typed fields, leading to extra type parameters.
unsafe_grid_getindex(g::IndexSubGrid, idx) = unsafe_grid_getindex(g.supergrid, g.subindices[idx])

function mask(g::IndexSubGrid)
    mask = falses(size(supergrid(g)))
	for i in g.subindices
		mask[i] = true
	end
    mask
end


covering(g::IndexSubGrid) = g.domain

# Check whether element grid[i] (of the underlying grid) is in the indexed subgrid.
issubindex(i, g::IndexSubGrid) = in(i, subindices(g))

getindex(grid::AbstractGrid, i::AbstractArray{Int}) = IndexSubGrid(grid, i)
