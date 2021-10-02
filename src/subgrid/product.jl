const ProductSubGrid = FlatProductGrid{N,T,<:NTuple{N,GRID}} where {N,T} where {GRID<:SubGrid}
# function tensorproductbitarray(vectors::Union{BitVector,Vector{Bool}}...)
# 	N = length(vectors)
# 	R = falses(map(length, vectors))
# 	for i in CartesianIndices(size(R))
# 		R[i] = reduce(&,map(k->vectors[k][i.I[k]],1:length(vectors)))
# 	end
# 	R
# end

function tensorproductbitarray(vectors::Vararg{<:Union{BitVector,Vector{Bool}}, N}) where N
	R = falses(map(length, vectors))
	@inbounds for i in CartesianIndices(size(R))
		t = true
		for k in 1:N
			if !vectors[k][i.I[k]]
				t = false
				break
			end
		end
		if t
			R[i] = t
		end
	end
	R
end

mask(grid::ProductSubGrid) = tensorproductbitarray(map(mask, components(grid))...)
subindices(grid::ProductSubGrid) = findall(mask(grid))
supergrid(grid::ProductSubGrid) = ProductGrid(map(supergrid, components(grid))...)
issubindex(i, g::ProductSubGrid) = all(map(issubindex, i, components(g)))
issubindex(i::CartesianIndex, g::ProductSubGrid) = issubindex(i.I, g)
