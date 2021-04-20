
abstract type AbstractMappedGrid{T,N} <: SimpleLazyGrid{T,N}
end

"""
A MappedGrid consists of a grid and a map. Each grid point of the mapped grid
is the map of the corresponding point of the underlying grid.
"""
struct MappedGrid{G,M,T,N} <: AbstractMappedGrid{T,N}
	grid	::	G
	map		::	M

	MappedGrid{G,M,T,N}(grid::AbstractGrid{T,N}, map) where {G,M,T,N} = new(grid, map)
end

const MappedGrid1d{G,M,T<:Number,N} = MappedGrid{G,M,T,N}

MappedGrid(grid::AbstractGrid{T,N}, map::AbstractMap) where {T,N} =
	MappedGrid{typeof(grid),typeof(map),T,N}(grid, map)

name(grid::MappedGrid) = "Mapped grid"

forward_map(g::MappedGrid) = g.map

mapped_grid(grid::AbstractGrid, map::AbstractMap) = MappedGrid(grid, map)

# avoid multiple mappings
mapped_grid(g::MappedGrid, map::AbstractMap) = MappedGrid(supergrid(g), map∘forward_map(g))

# Convenience function, similar to apply_map for Dictionary's
apply_map(grid::AbstractGrid, map::AbstractMap) = mapped_grid(grid, map)

for op in (:length, :size, :eachindex, :indextype, :isperiodic)
	@eval $op(g::MappedGrid) = $op(supergrid(g))
end

for op in (:minimum, :maximum, :covering)
	@eval $op(g::MappedGrid1d) = forward_map(g).($op(supergrid(g)))
end

resize(g::MappedGrid, n::Int) = apply_map(resize(supergrid(g), n), forward_map(g))

unsafe_grid_getindex(g::MappedGrid, idx...) = g.map(g.grid[idx...])

function rescale(g::AbstractGrid1d, a, b)
	m = mapto(covering(g), a..b)
	mapped_grid(g, m)
end


# Preserve tensor product structure
function rescale(g::ProductGrid, a::SVector{N}, b::SVector{N}) where {N}
	scaled_grids = [ rescale(component(g, i), a[i], b[i]) for i in 1:N]
	ProductGrid(scaled_grids...)
end
