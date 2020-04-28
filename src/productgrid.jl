
"""
A `ProductGrid` represents the cartesian product of other grids.

`struct ProductGrid{TG,T,N} <: AbstractGrid{T,N}`

Parameters:
- TG is a tuple of (grid) types
- T is the element type of the grid
- N is the dimension of the grid layout
"""
struct ProductGrid{TG,T,N} <: AbstractGrid{T,N}
	grids	::	TG
end

# Generic functions for composite types:
elements(grid::ProductGrid) = grid.grids
element(grid::ProductGrid, j::Int) = grid.grids[j]
element(grid::ProductGrid, range::AbstractRange) = cartesianproduct(grid.grids[range]...)
iscomposite(::ProductGrid) = true

function ProductGrid(grids...)
	T = mapreduce(numtype, promote_type, grids)
	N = sum(map(dimension, grids))
	ProductGrid{typeof(grids),SVector{N,T},length(grids)}(grids)
end

size(g::ProductGrid) = map(length, g.grids)
size(g::ProductGrid, j::Int) = length(g.grids[j])

gridsupport(g::ProductGrid) = cartesianproduct(map(gridsupport, elements(g))...)
isperiodic(g::ProductGrid) = reduce(&, map(isperiodic, elements(g)))

getindex(g::ProductGrid{TG,T,N}, I::Vararg{Int,N}) where {TG,T,N} =
	convert(T, map(getindex, g.grids, I))

similargrid(grid::ProductGrid, ::Type{T}, dims...) where T = error()#ProductGrid([similargrid(g, eltype(T), dims[i]) for (i,g) in enumerate(elements(grid))]...)

toiterator(::Type{ProductGrid}, a::AbstractGrid) = (a,)
toiterator(::Type{ProductGrid}, a::ProductGrid) = Iterators.flatten(ntuple(k->toiterator(ProductGrid,element(a,k)),Val(dimension(a))))
toiterator(::Type{ProductGrid}, elements::Vararg{AbstractGrid,N}) where N = Iterators.flatten(ntuple(k->toiterator(ProductGrid,elements[k]),Val(N)))
flatten(::Type{ProductGrid}, elements::Vararg{AbstractGrid,N}) where N = flatten(ProductGrid, toiterator(ProductGrid, elements...))
function flatten(::Type{ProductGrid}, a::Iterators.Flatten)
	collect(a)
end

for (BaseType,TPType) in [ (:AbstractGrid, :ProductGrid)]
    # Override × for grids
    @eval cross(args::$BaseType...) = cartesianproduct(args...)
	@eval ^(arg::$BaseType, n::Int) = cartesianproduct(arg, n)
    # In order to avoid strange nested structures, we flatten the arguments
    @eval cartesianproduct(args::$BaseType...) = $TPType(flatten($TPType, args...)...)
    @eval cartesianproduct(arg::$BaseType, n::Int) = cartesianproduct(ntuple(k->arg,Val(n)))
    # Disallow cartesian products with just one argument
    @eval cartesianproduct(arg::$BaseType) = arg
end
