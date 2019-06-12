
@recipe function f(grid::AbstractGrid{T,N}, vals::AbstractMatrix) where {T,N}
    seriestype --> :surface
    size --> (500,400)
    xrange = LinRange(endpoints(support(element(grid,1)))..., size(grid,1))
    yrange = LinRange(endpoints(support(element(grid,2)))..., size(grid,2))
    xrange, yrange, vals'
end

@recipe function f(grid::AbstractGrid{T,N}, vals::AbstractVector) where {T,N}
    seriestype --> :surface
    legend=false
    size --> (500,400)
    xrange = [x[1] for x in grid]
    yrange = [x[2] for x in grid]
    xrange, yrange, vals
end

# Plot an Nd grid
@recipe function f(grid::AbstractGrid)
    seriestype --> :scatter
    size --> (500,400)
    legend --> false
    broadcast(x->tuple(x...),collect(grid)[1:end])
end

# Plot a 1D grid
@recipe function f(grid::AbstractGrid1d)
    seriestype --> :scatter
    yticks --> []
    ylims --> [-1 1]
    size --> (500,200)
    legend --> false
    collect(grid), zeros(size(grid))
end

# # Plot a matrix of values on a 2D equispaced grid
# @recipe function f(grid::AbstractGrid2d, vals)
#     seriestype --> :surface
#     size --> (500,400)
#     xrange = linspace(leftendpoint(grid)[1],rightendpoint(grid)[1],size(grid,1))
#     yrange = linspace(leftendpoint(grid)[2],rightendpoint(grid)[2],size(grid,2))
#     xrange, yrange, vals'
# end
