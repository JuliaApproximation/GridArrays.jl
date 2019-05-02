# from grid/grid.jl
module Grids
using DomainSets, StaticArrays, RecipesBase, Test

using DomainSets: endpoints

import Base: *, size, length, @propagate_inbounds, step, ndims, unsafe_getindex,
    checkbounds, IndexStyle, ==, ≈, getindex, eachindex, convert, in, ^
import Base.Broadcast: broadcast

import DomainSets: cartesianproduct, element, elements, numelements, ×, cross, minimum, maximum


export AbstractGrid, AbstractGrid1d, AbstractGrid3d,
        AbstractEquispacedGrid, EquispacedGrid, PeriodicEquispacedGrid,
        FourierGrid, MidpointEquispacedGrid, RandomEquispacedGrid,
        AbstractIntervalGrid, eachelement, ScatteredGrid, ×, cartesianproduct,
        TensorSubGrid, instantiate, support, float_type, isperiodic, MaskedGrid,
        boundary, subgrid, mask, randomgrid, boundingbox
export ChebyshevNodes, ChebyshevGrid, ChebyshevPoints, ChebyshevExtremae, ×
export Point

# from grid/productgrid.jl
export ProductGrid

# from grid/subgrid.jl
export subindices, supergrid, issubindex, similar_subgrid

# from grid/mappedgrid.jl
export MappedGrid, mapped_grid, apply_map


"Assign a floating point type to a domain element type T."
float_type(::Type{T}) where {T <: Real} = T
float_type(::Type{Complex{T}}) where {T <: Real} = Complex{T}
float_type(::Type{SVector{N,T}}) where {N,T} = T
float_type(::Type{NTuple{N,T}}) where {N,T} = T

# Fallback: we return Float64
float_type(::Type{T}) where {T} = Float64


include("domains/extensions.jl")

include("grid.jl")
include("productgrid.jl")
include("intervalgrids.jl")
include("mappedgrid.jl")
include("scattered_grid.jl")
include("randomgrid.jl")


include("subgrid/AbstractSubGrids.jl")


include("recipes.jl")


export test_generic_grid, test_interval_grid
export grid_iterator1, grid_iterator2
include("test/test_grids.jl")


≈(d1::DomainSets.Interval, d2::DomainSets.Interval) =
    1≈1+abs(DomainSets.leftendpoint(d1)-DomainSets.leftendpoint(d2))+abs(DomainSets.rightendpoint(d1)-DomainSets.rightendpoint(d2))


end # module
