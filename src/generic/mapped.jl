
"An `AbstractMappedGrid` represents the lazy application of a map to a grid."
abstract type AbstractMappedGrid{T,N} <: SimpleLazyGrid{T,N} end

const AbstractMappedGrid1d{T<:Number,N} = AbstractMappedGrid{T,N}

for op in (:length, :size, :eachindex, :indextype, :isperiodic)
	@eval $op(g::AbstractMappedGrid) = $op(supergrid(g))
end

for op in (:minimum, :maximum, :covering)
	@eval $op(g::AbstractMappedGrid1d) = forward_map(g).($op(supergrid(g)))
end

resize(g::AbstractMappedGrid, n::Int) = map_grid(resize(supergrid(g), n), forward_map(g))

function rescale(g::AbstractGrid1d, a, b)
	m = mapto(covering(g), a..b)
	map_grid(g, m)
end

map_grid(grid::AbstractGrid, map) = map_grid1(grid, map)
map_grid1(grid::AbstractGrid, map) = map_grid2(grid, map)
map_grid2(grid, map) = MappedGrid(grid, map)

# some simplifications
map_grid1(g::AbstractMappedGrid, map) = map_grid(supergrid(g), map∘forward_map(g))
map_grid2(grid, map::IdentityMap) = grid

# Convenience function, similar to apply_map for Dictionary's
apply_map(grid::AbstractGrid, map::AbstractMap) = map_grid(grid, map)

forward_map(g::AbstractGrid, x...) = forward_map(g)(x...)
inverse_map(g::AbstractGrid, x...) = inverse_map(g)(x...)

unsafe_grid_getindex(g::AbstractMappedGrid, I...) = forward_map(g, supergrid(g, I...))


# Preserve tensor product structure
function rescale(g::ProductGrid, a::SVector{N}, b::SVector{N}) where {N}
	scaled_grids = [ rescale(component(g, i), a[i], b[i]) for i in 1:N]
	ProductGrid(scaled_grids...)
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

MappedGrid(grid::AbstractGrid{T,N}, map) where {T,N} =
	MappedGrid{typeof(grid),typeof(map),T,N}(grid, map)

name(grid::MappedGrid) = "Mapped grid"

forward_map(g::MappedGrid) = g.map
inverse_map(g::MappedGrid) = inverse(forward_map(g))
inverse_map(g::MappedGrid, x) = inverse(forward_map(g), x)
