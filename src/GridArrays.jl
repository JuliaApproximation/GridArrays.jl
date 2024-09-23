module GridArrays

using CompositeTypes,
        CompositeTypes.Display,
        DomainSets,
        FastGaussQuadrature,
        FillArrays,
        GaussQuadrature, 
        StaticArrays,
        Test

using DomainSets: endpoints

## List of imports

import Base:
        size,
        length,
        @propagate_inbounds,
        range,
        step,
        ndims,
        checkbounds,
        IndexStyle,
        ==, ≈,
        getindex, eachindex,
        convert, in, ^, string, axes,
        show

import CompositeTypes:
        iscomposite, component, components, ncomponents

import DomainSets:
        factors,
        ×, cross,
        dimension, prectype, numtype,
        minimum, maximum,
        boundary,
        mapto_canonical, mapfrom_canonical,
        forward_map, inverse_map


## List of exports

# from generic/grid.jl
export AbstractGrid, AbstractGrid1d,
        numtype, prectype,
        dimension,
        resize,
        covering,
        isperiodic,
        boundary,
        canonicalgrid

# from generic/product.jl
export ProductGrid, productgrid,
        VcatGrid, FlatProductGrid

# from subgrid/subgrid.jl
export subindices, supergrid, issubindex, similar_subgrid

# from generic/mapped.jl
export MappedGrid, map_grid

# from subgrid/subgrid.jl
export subgrid, mask, ProductSubGrid

# from domains/interval.jl
export AbstractIntervalGrid,
        AbstractEquispacedGrid,
        EquispacedGrid, PeriodicEquispacedGrid, MidpointEquispacedGrid,
        UnitEquispacedGrid, UnitPeriodicEquispacedGrid, UnitMidpointEquispacedGrid,
        FourierGrid,
        step

# from domains/scattered.jl
export ScatteredGrid, randomgrid

# from domains/discretize.jl
export discretize_togrid

# from applications/gauss.jl
export ChebyshevTNodes, ChebyshevUNodes, ChebyshevNodes,
        ChebyshevExtremae,
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
include("domains/discretize.jl")

include("subgrid/subgrid.jl")

include("applications/gauss.jl")

# We define some testing routines inside the package, so that they
# can also be used in the tests of other packages that extend grids
include("test/Test.jl")


export coverdomain
coverdomain(g) = covering(g)


end # module
