module GridArrays

using DomainSets, StaticArrays, RecipesBase, Test, FastGaussQuadrature,
        GaussQuadrature, FillArrays

using DomainSets: endpoints

import Base: *, size, length, @propagate_inbounds, step, ndims,
    checkbounds, IndexStyle, ==, ≈, getindex, eachindex, convert, in, ^, string, axes, convert
import Base.Broadcast: broadcast

import DomainSets: cartesianproduct, iscomposite, element, elements, numelements, ×, cross,
        minimum, maximum,
        dimension, prectype, numtype


export AbstractGrid, AbstractGrid1d, AbstractGrid3d,
        AbstractEquispacedGrid, EquispacedGrid, PeriodicEquispacedGrid,
        FourierGrid, MidpointEquispacedGrid, RandomEquispacedGrid,
        AbstractIntervalGrid, eachelement, ScatteredGrid, cartesianproduct,
        TensorSubGrid, covering, isperiodic,
        boundary, subgrid, mask, randomgrid, boundingbox, TensorSubGrid,
        dimension, prectype, ChebyshevTNodes, ChebyshevUNodes,
        LaguerreNodes, HermiteNodes, LegendreNodes, JacobiNodes
export ChebyshevNodes, ChebyshevGrid, ChebyshevPoints, ChebyshevExtremae

# from productgrid.jl
export ProductGrid

# from subgrid/subgrid.jl
export subindices, supergrid, issubindex, similar_subgrid

# from mappedgrid.jl
export MappedGrid, mapped_grid, apply_map

import Base: isapprox

# TODO move to domainsets
isapprox(d1::DomainSets.ProductDomain,d2::DomainSets.ProductDomain) =
    reduce(&, map(isapprox,DomainSets.elements(d1), DomainSets.elements(d2)))


include("util/PeriodicCartesianIndices.jl")
using ..PeriodicIndexing


include("grid.jl")
include("productgrid.jl")
include("intervalgrids.jl")
include("mappedgrid.jl")
include("scattered_grid.jl")

include("domains/boundingbox.jl")
include("domains/broadcast.jl")

include("subgrid/AbstractSubGrids.jl")

include("applications/gauss.jl")
include("applications/randomgrid.jl")

include("recipes.jl")

# We define some testing routines inside the package, so that they
# can also be used in the tests of other packages that extend grids
include("test/Test.jl")


≈(d1::DomainSets.Interval, d2::DomainSets.Interval) =
    1≈1+abs(DomainSets.infimum(d1)-DomainSets.infimum(d2))+abs(DomainSets.supremum(d1)-DomainSets.supremum(d2))

export coverdomain
coverdomain(g) = covering(g)



end # module
