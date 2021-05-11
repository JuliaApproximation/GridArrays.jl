
"An `AbstractMappedGrid` represents the lazy application of a map to a grid."
abstract type AbstractMappedGrid{T,N} <: SimpleLazyGrid{T,N} end

const AbstractMappedGrid1d{T<:Number,N} = AbstractMappedGrid{T,N}

for op in (:length, :size, :eachindex, :indextype, :isperiodic)
	@eval $op(g::AbstractMappedGrid) = $op(supergrid(g))
end

for op in (:minimum, :maximum, :covering)
	@eval $op(g::AbstractMappedGrid1d) = forward_map(g).($op(supergrid(g)))
end

resize(g::AbstractMappedGrid, n::Int) = map_grid(forward_map(g), resize(supergrid(g), n))

function rescale(g::AbstractGrid1d, a, b)
	m = mapto(covering(g), a..b)
	map_grid(m, g)
end

(→)(g::AbstractGrid1d, d::AbstractInterval) = rescale(g, infimum(d), supremum(d))

map_grid(grid::AbstractGrid, map) = (@warn "Interchange arguments of map_grid here"; map_grid(map, grid))
map_grid(map, grid::AbstractGrid) = map_grid1(map, grid)
map_grid1(map, grid::AbstractGrid) = map_grid2(map, grid)
map_grid2(map, grid) = MappedGrid(map, grid)

# some simplifications
map_grid1(map, g::AbstractMappedGrid) = map_grid(map∘forward_map(g), supergrid(g))
map_grid2(map::IdentityMap, grid) = grid

# Convenience function, similar to apply_map for Dictionary's
apply_map(grid::AbstractGrid, map::AbstractMap) = map_grid(map, grid)

forward_map(g::AbstractGrid, x...) = forward_map(g)(x...)
inverse_map(g::AbstractGrid, x...) = inverse_map(g)(x...)

unsafe_grid_getindex(g::AbstractMappedGrid, I...) = forward_map(g, supergrid(g, I...))


# Preserve tensor product structure
function rescale(g::ProductGrid, a::SVector{N}, b::SVector{N}) where {N}
	scaled_grids = [ rescale(component(g, i), a[i], b[i]) for i in 1:N]
	ProductGrid(scaled_grids...)
end


# Base.show(io::IO, mime::MIME"text/plain", g::AbstractMappedGrid) =
# 	composite_show(io, mime, g)
# Display.displaystencil(g::AbstractMappedGrid) =
#     DomainSets.map_stencil_broadcast(forward_map(g), supergrid(g))
# Display.object_parentheses(g::AbstractMappedGrid) =
#     Display.object_parentheses(forward_map(g))
# Display.stencil_parentheses(g::AbstractMappedGrid) =
#     Display.stencil_parentheses(forward_map(g))


"""
A MappedGrid consists of a grid and a map. Each grid point of the mapped grid
is the map of the corresponding point of the underlying grid.
"""
struct MappedGrid{T,N,M,G} <: AbstractMappedGrid{T,N}
	map		::	M
	grid	::	G
end

const MappedGrid1d{T<:Number,N,M,G} = MappedGrid{T,N,M,G}

MappedGrid(map, grid::AbstractGrid{T,N}) where {T,N} =
	MappedGrid{T,N}(map, grid)
MappedGrid(map::Map{S}, grid::AbstractGrid{T,N}) where {S,T,N} =
	MappedGrid{codomaintype(map),N}(map, grid)
MappedGrid{T,N}(map, grid) where {T,N} =
	MappedGrid{T,N,typeof(map),typeof(grid)}(map, grid)

MappedGrid(grid::AbstractGrid, map) = MappedGrid(map, grid)

forward_map(g::MappedGrid) = g.map
inverse_map(g::MappedGrid) = inverse(forward_map(g))
inverse_map(g::MappedGrid, x) = inverse(forward_map(g), x)
