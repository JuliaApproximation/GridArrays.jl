module GridArraysRecipesBaseExt

using GridArrays, RecipesBase

using GridArrays:
    component,
    endpoints

@recipe function f(grid::AbstractGrid{T,N}, vals::AbstractMatrix) where {T,N}
    seriestype --> :surface
    size --> (500,400)
    xrange = LinRange(endpoints(covering(component(grid,1)))..., size(grid,1))
    yrange = LinRange(endpoints(covering(component(grid,2)))..., size(grid,2))
    xrange, yrange, vals'
end

@recipe function f(grid::AbstractGrid{<:AbstractVector}, vals::AbstractVector)
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

# Disregarded as this recipe does not involve a type specific to GridArrays
# @recipe function f(a::Vector{SArray{Tuple{2},T,1,2}}) where T
#     A = zeros(T,length(a),2)
#     for (l,e) in enumerate(a)
#             A[l,:] .= e
#     end
#     A[:,1],A[:,2]
# end

end
