module GridArrays

using DomainSets, StaticArrays, RecipesBase, Test, FastGaussQuadrature,
        GaussQuadrature, FillArrays

using DomainSets: endpoints

@deprecate AbstractSubGrid SubGrid
@deprecate MaskedSubGrid MaskedGrid
@deprecate mapping = forward_map

## List of imports

import Base: size, length, @propagate_inbounds, step, ndims,
    checkbounds, IndexStyle, ==, ≈, getindex, eachindex,
    convert, in, ^, string, axes

import DomainSets:
        iscomposite, component, components, ncomponents,
        ×, cross,
        minimum, maximum,
        boundary,
        dimension, prectype, numtype,
        tocanonical, fromcanonical


## List of exports

# from util/common.jl
export numtype, prectype

# from generic/grid.jl
export AbstractGrid, AbstractGrid1d, AbstractGrid3d,
        dimension,
        covering,
        isperiodic,
        boundary

# from generic/product.jl
export ProductGrid, productgrid,
        VcatGrid, FlatProductGrid

# from subgrid/subgrid.jl
export subindices, supergrid, issubindex, similar_subgrid

# from generic/mapped.jl
export MappedGrid, mapped_grid, apply_map

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
