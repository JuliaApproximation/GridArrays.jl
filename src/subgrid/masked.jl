
# See also the file grid/subgrid.jl in BasisFunctions for the definition of
# SubGrid and IndexSubGrid.


"""
    struct MaskedGrid{T,GRID,MASK,I,D} <: SubGrid{T,1}

A MaskedGrid is a subgrid of another grid that is defined by a mask.
The mask is true or false for each point in the supergrid. The set of points
for which it is true make up the MaskedGrid.
"""
struct MaskedGrid{T,GRID,MASK,I,D} <: SubGrid{T,1}
    supergrid   ::	GRID
    mask	    ::	MASK
    indices     ::  Vector{I}
    M           ::	Int				# Total number of points in the mask
    domain      ::  D

    MaskedGrid{T,GRID,MASK,I,D}(supergrid::AbstractArray{T}, mask, indices, domain) where {T,GRID,MASK,I,D} =
        new(supergrid, mask, indices, sum(mask), domain)
end
# TODO: In MaskedGrid, perhaps we should not be storing pointers to the points of the underlying grid, but
# rather the points themselves. In that case we wouldn't need to specialize on the type of grid (parameter G can go).


MaskedGrid(supergrid::AbstractArray{T}, mask, indices, domain) where T =
	MaskedGrid{T}(supergrid, mask, indices, domain)

function MaskedGrid{T}(supergrid, mask, indices, domain) where T
	@assert size(supergrid) == size(mask)
	MaskedGrid{T,typeof(supergrid),typeof(mask),eltype(indices),typeof(domain)}(supergrid, mask, indices, domain)
end

# These are for the assignment to indices in the function below.
convert(::Type{NTuple{N,Int}},i::CartesianIndex{N}) where {N} = ntuple(k->i[k],N)

MaskedGrid(supergrid::AbstractArray, domain::Domain) =
	MaskedGrid(supergrid, in.(supergrid, Ref(domain)), domain)

# MaskedGrid(maskedgrid::MaskedGrid, domain::Domain) =
#     MaskedGrid(supergrid(maskedgrid), mask(maskedgrid) .& in.(supergrid(maskedgrid), domain))

MaskedGrid(supergrid::AbstractArray, mask, domain) =
    MaskedGrid(supergrid, mask, subindices(supergrid, mask), domain)

function subindices(supergrid, mask::BitArray)
    I = eltype(eachindex(supergrid))
    indices = Array{I}(undef, sum(mask))
    i = 1
    for m in eachindex(supergrid)
       if mask[m]
           indices[i] = m
           i += 1
       end
    end
    indices
end

size(g::MaskedGrid) = (g.M,)

mask(g::MaskedGrid) = g.mask

covering(g::MaskedGrid) = g.domain

subindices(g::MaskedGrid) = g.indices

similar_subgrid(g::MaskedGrid, g2::AbstractArray) = MaskedGrid(g2, g.mask, g.indices)


# Check whether element grid[i] (of the underlying grid) is in the masked grid.
issubindex(i, g::MaskedGrid) = g.mask[i]

unsafe_grid_getindex(g::MaskedGrid, idx::Int) = @inbounds getindex(g.supergrid, g.indices[idx])

getindex(g::AbstractGrid, idx::BitArray) = MaskedGrid(g, idx, covering(g))

function subgrid(grid::MaskedGrid, domain::Domain)
    submask = in.(supergrid(grid), Ref(domain))
    MaskedGrid(supergrid(grid), submask .& mask(grid), domain)
end
