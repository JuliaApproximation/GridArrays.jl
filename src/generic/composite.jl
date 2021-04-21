
"Supertype of lazy grids, defined in terms of other grids."
abstract type LazyGrid{T,N} <: AbstractGrid{T,N}
end

"Grid defined in terms of a single `supergrid`."
abstract type SimpleLazyGrid{T,N} <: AbstractGrid{T,N}
end

supergrid(g::SimpleLazyGrid) = g.grid
supergrid(g::SimpleLazyGrid, I...) = supergrid(g)[I...]

"Grid defined in terms of multiple grids."
abstract type CompositeGrid{T,N} <: LazyGrid{T,N}
end

components(grid::CompositeGrid) = grid.grids
