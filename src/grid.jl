
"""
    abstract type AbstractGrid{T,N} <: AbstractArray{T,N}

Grids are arrays of points.
"""
abstract type AbstractGrid{T,N} <: AbstractArray{T,N}
end

const AbstractGrid1d{T <: Real} = AbstractGrid{T,1}

prectype(::Type{G}) where {G<:AbstractGrid} = prectype(eltype(G))
numtype(::Type{G}) where {G<:AbstractGrid} = numtype(eltype(G))

dimension(::AbstractGrid{<:Number}) = 1
dimension(::AbstractGrid{<:SVector{N,T}}) where {N,T} = N
dimension(::AbstractGrid{<:NTuple{N,Any}}) where {N} = N

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

@deprecate support(grid::AbstractGrid) coverdomain(grid) false

export resize
"""
    resize(grid::AbstractGrid, dims...)

Create a grid of same structure (such as product structure)
but with different points in the different dimensions.
"""
resize(grid::AbstractGrid{T}, dims...) where {T} = similargrid(grid, T, dims...)

"""
    hasextension(grid::AbstractGrid)

Is it possible to use the `resize` function?
See also [`resize`](@ref)
"""
hasextension(::AbstractGrid) = false
