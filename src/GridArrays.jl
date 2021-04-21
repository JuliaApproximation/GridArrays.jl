module GridArrays

using CompositeTypes, DomainSets, StaticArrays, RecipesBase, Test, FastGaussQuadrature,
        GaussQuadrature, FillArrays
using CompositeTypes.Display

using DomainSets: endpoints

@deprecate AbstractSubGrid SubGrid
@deprecate MaskedSubGrid MaskedGrid
@deprecate mapping forward_map
@deprecate mapped_grid map_grid

## List of imports

import Base: size, length, @propagate_inbounds, step, ndims,
    checkbounds, IndexStyle, ==, ≈, getindex, eachindex,
    convert, in, ^, string, axes

import DomainSets:
        iscomposite, component, components, ncomponents,
        ×, cross,
        dimension, prectype, numtype,
        minimum, maximum,
        boundary,
        tocanonical, fromcanonical,
        forward_map, inverse_map


## List of exports

# from generic/grid.jl
export AbstractGrid, AbstractGrid1d, AbstractGrid3d,
        numtype, prectype,
        dimension,
        resize,
        covering,
        isperiodic,
        boundary

# from generic/product.jl
export ProductGrid, productgrid,
        VcatGrid, FlatProductGrid

# from subgrid/subgrid.jl
export subindices, supergrid, issubindex, similar_subgrid

# from generic/mapped.jl
export MappedGrid, map_grid, apply_map

# from subgrid/subgrid.jl
export subgrid, mask, TensorSubGrid

# from domains/interval.jl
export AbstractIntervalGrid,
        AbstractEquispacedGrid, EquispacedGrid, PeriodicEquispacedGrid,
        FourierGrid, MidpointEquispacedGrid, RandomEquispacedGrid

# from domains/scattered.jl
export ScatteredGrid, randomgrid

# from applications/gauss.jl
export ChebyshevTNodes, ChebyshevUNodes, ChebyshevNodes,
        ChebyshevGrid, ChebyshevPoints, ChebyshevExtremae,
        LaguerreNodes, HermiteNodes, LegendreNodes, JacobiNodes


include("util/PeriodicCartesianIndices.jl")
using ..PeriodicIndexing

include("util/common.jl")

include("generic/grid.jl")
include("generic/composite.jl")
include("generic/product.jl")
include("generic/mapped.jl")

include("domains/interval.jl")
include("domains/scattered.jl")

include("subgrid/subgrid.jl")

include("applications/gauss.jl")

include("util/recipes.jl")

# We define some testing routines inside the package, so that they
# can also be used in the tests of other packages that extend grids
include("test/Test.jl")


export coverdomain
coverdomain(g) = covering(g)


end # module
