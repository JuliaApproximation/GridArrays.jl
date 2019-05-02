"""
    abstract type AbstractSubGrid{T,N} <: AbstractGrid{T,N} end

A subgrid of an underlying grid.
"""
abstract type AbstractSubGrid{T,N} <: AbstractGrid{T,N} end

"""
    supergrid(g::AbstractSubGrid)

The underlying grid of the subgrid.
"""
supergrid(g::AbstractSubGrid) = g.supergrid

include("indexsubgrid.jl")
include("maskedsubgrid.jl")
