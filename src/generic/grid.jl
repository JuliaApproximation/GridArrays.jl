
"""
    abstract type AbstractGrid{T,N} <: AbstractArray{T,N}

Grids are arrays of points.
"""
abstract type AbstractGrid{T,N} <: AbstractArray{T,N}
end

const AbstractGridVector{T} = AbstractGrid{T,1}
const AbstractGrid1d{T <: Number} = AbstractGrid{T,1}

const GridLike{T,N} = Union{AbstractGrid{T,N},AbstractArray{T,N}}
const Grid1dLike{T<:Number,N} = Union{AbstractGridVector{T},AbstractVector{T}}


prectype(::Type{G}) where {G<:AbstractGrid} = prectype(eltype(G))
numtype(::Type{G}) where {G<:AbstractGrid} = numtype(eltype(G))
dimension(::AbstractGrid{T}) where {T} = DomainSets.euclideandimension(T)

# Do a bounds check and invoke unsafe_grid_getindex.
# Concrete subtypes can specialize unsafe_grid_getindex without checking bounds.
@propagate_inbounds function getindex(grid::AbstractGrid{T,1}, i::Int) where {T}
    checkbounds(grid, i)
    unsafe_grid_getindex(grid, i)
end

@propagate_inbounds function getindex(grid::AbstractGrid{T,N}, i::Vararg{Int,N}) where {T,N}
    checkbounds(grid, i...)
    unsafe_grid_getindex(grid, i...)
end

convert(::Type{AbstractGrid{T}}, grid::AbstractGrid{T,N}) where {T,N} = grid
convert(::Type{AbstractGrid{T}}, grid::AbstractGrid{S,N}) where {S,T,N} = similargrid(grid, S, size(grid))

similargrid(grid::AbstractGrid1d, ::Type{T}, dims::Tuple{Int}) where {T} = similargrid(grid, T, dims[1])


"""
    resize(grid::AbstractGrid, dims...)

Create a grid of the same structure (such as product structure),
but with different points and different dimensions.
"""
resize(grid::AbstractGrid{T}, dims...) where {T} = similargrid(grid, T, dims...)

"""
    hasextension(grid::AbstractGrid)

Is it possible to use the `resize` function?
See also [`resize`](@ref)
"""
hasextension(::AbstractGrid) = false

"Return the canonical grid, if any, of the family of grids."
canonicalgrid(g::AbstractGrid) = g

Display.displaysymbol(g::AbstractGrid) = 'g'
