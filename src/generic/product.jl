
"A `ProductGrid` represents the tensor product of other grids."
abstract type ProductGrid{T,N} <: CompositeGrid{T,N}
end

size(g::ProductGrid) = map(length, g.grids)
size(g::ProductGrid, j::Int) = length(g.grids[j])

ProductGrid(grids...) = VcatGrid(grids...)
ProductGrid(grids::Grid1dLike...) = FlatProductGrid(grids...)

ProductGrid{T}(grids...) where {N,S,T <: SVector{N,S}} = VcatGrid{N,S}(grids...)
ProductGrid{T}(grids::Grid1dLike...) where {N,S,T <: SVector{N,S}} =
	FlatProductGrid{N,S}(grids...)

similargrid(grid::ProductGrid, ::Type{T}, dims...) where T = error()#ProductGrid([similargrid(g, eltype(T), dims[i]) for (i,g) in enumerate(components(grid))]...)
resize(grid::ProductGrid, n) = ProductGrid(map(resize, components(grid), n)...)

covering(g::ProductGrid) = productdomain(map(covering, components(g))...)

isperiodic(g::ProductGrid) = mapreduce(isperiodic, &, components(g))

factors(g::ProductGrid) = components(g)

unsafe_grid_getindex(g::ProductGrid{T,N}, I::Vararg{Int,N}) where {T,N} =
	convert(T, map(getindex, components(g), I))

productgrid() = ()
productgrid(d) = d
productgrid(d1, d2, d3...) = productgrid(productgrid(d1, d2), d3...)

productgrid(d1, d2) = productgrid1(d1, d2)
productgrid1(d1, d2) = productgrid2(d1, d2)
productgrid2(d1, d2) = ProductGrid(d1, d2)

productgrid(d1::ProductGrid, d2::ProductGrid) =
	ProductGrid(components(d1)..., components(d2)...)
productgrid1(d1::ProductGrid, d2) = ProductGrid(components(d1)..., d2)
productgrid2(d1, d2::ProductGrid) = ProductGrid(d1, components(d2)...)

cross(grids::AbstractGrid...) = productgrid(grids...)
^(grid::AbstractGrid, n::Int) = productgrid(ntuple(i->grid, n)...)

canonicalgrid(g::ProductGrid) = ProductGrid(map(canonicalgrid, components(g)))

mapto_canonical(g::ProductGrid) = ProductMap(map(mapto_canonical, components(g)))
mapfrom_canonical(g::ProductGrid) = ProductMap(map(mapfrom_canonical, components(g)))

# Pretty printing
Display.combinationsymbol(d::ProductGrid) = Display.Times()
Display.displaystencil(d::ProductGrid) = composite_displaystencil(d)
Base.show(io::IO, mime::MIME"text/plain", d::ProductGrid) = composite_show(io, mime, d)
Base.show(io::IO, d::ProductGrid) = composite_show_compact(io, d)


"A `FlatProductGrid` is a product grid of `N` 1-D grids."
struct FlatProductGrid{N,T,GG} <: ProductGrid{SVector{N,T},N}
	grids	::	GG
end

FlatProductGrid(grids::Vararg{Any,N}) where {N} =
	FlatProductGrid{N}(grids...)
function FlatProductGrid{N}(grids::Vararg{Any,N}) where {N}
	T = mapreduce(numtype, promote_type, grids)
	FlatProductGrid{N,T}(grids...)
end
FlatProductGrid{N,T}(grids::Vararg{Any,N}) where {N,T} =
	FlatProductGrid{N,T,typeof(grids)}(grids)


"""
A `VcatGrid` is a product grid that concatenates the elements of all composing
grids into a flat vector, like `vcat` does.
"""
struct VcatGrid{N,T,M,DIMS,GG} <: ProductGrid{SVector{N,T},M}
	grids	::	GG
end

function VcatGrid(grids::Vararg{AbstractGrid{<:Any,<:Any},M}) where {M}
	T = mapreduce(numtype, promote_type, grids)
	DIMS = map(dimension, grids)
	N = sum(DIMS)
	VcatGrid{N,T,M,DIMS}(grids...)
end
function VcatGrid{N,T}(grids::Vararg{AbstractGrid{<:Any,<:Any},M}) where {N,T,M}
	DIMS = map(dimension, grids)
	VcatGrid{N,T,M,DIMS}(grids...)
end

VcatGrid{N,T,M,DIMS}(grids...) where {N,T,M,DIMS} =
	VcatGrid{N,T,M,DIMS,typeof(grids)}(grids)

unsafe_grid_getindex(g::VcatGrid{N,T,M,DIMS}, I::Vararg{Int,M}) where {N,T,M,DIMS} =
	DomainSets.convert_tocartesian(map(getindex, components(g), I), Val(DIMS))
