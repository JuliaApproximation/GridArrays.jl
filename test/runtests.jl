
using GridArrays

using DomainSets, FunctionMaps, LinearAlgebra, Plots, StaticArrays
using Test

using DomainSets: ×

using GridArrays: MaskedGrid, IndexSubGrid, randompoint,
    component, components, productgrid

using GridArrays.Test: test_generic_grid, test_interval_grid,
    grid_iterator1, grid_iterator2


function delimit(s::AbstractString)
    println()
    println("## ",s)
end

types = (Float64,BigFloat)


include("test_grids.jl")
include("test_gauss.jl")
include("test_periodiccartesianindices.jl")
include("test_boundingbox.jl")
include("test_gauss.jl")
include("test_subgrids.jl")
include("test_randomgrid.jl")
include("test_plots.jl")
