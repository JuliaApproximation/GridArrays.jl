# Plot a matrix of values on a 2D equispaced grid
@recipe function f(grid::AbstractGrid{T,N}, vals) where {T,N}
    seriestype --> :surface
    size --> (500,400)
    xrange = LinRange(minimum(grid)[1],maximum(grid)[1],size(grid,1))
    yrange = LinRange(minimum(grid)[2],maximum(grid)[2],size(grid,2))
    xrange, yrange, vals'
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
