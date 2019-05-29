var documenterSearchIndex = {"docs":
[{"location":"#Grids.jl-Documentation-1","page":"Home","title":"Grids.jl Documentation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"For installation instructions, see Installation.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"For a  full description of the functionality use the manual:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\"man/Grids.md\"]","category":"page"},{"location":"#Installation-1","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Grids.jl is not added to the Julia General registry.","category":"page"},{"location":"#(Recomanded)-1","page":"Home","title":"(Recomanded)","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"For Julia 1.1 or higher, you can add the FrameFun registry and than add Grids. From the Julia REPL, type ] to enter Pkg mode and run","category":"page"},{"location":"#","page":"Home","title":"Home","text":"pkg> add https://github.com/vincentcp/FrameFunRegistry\npkg> add BasisFunctions","category":"page"},{"location":"#Legacy-1","page":"Home","title":"Legacy","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"In Julia 1.0, the packages can be installed by cloning their git repository. From the Julia REPL, type ] to enter Pkg mode and run","category":"page"},{"location":"#","page":"Home","title":"Home","text":"pkg> add https://github.com/vincentcp/Grids.jl","category":"page"},{"location":"#","page":"Home","title":"Home","text":"or in a file you could use","category":"page"},{"location":"#","page":"Home","title":"Home","text":"using Pkg\npkg\"add https://github.com/vincentcp/Grids.jl\"","category":"page"},{"location":"man/Grids/#Grids-1","page":"Manual","title":"Grids","text":"","category":"section"},{"location":"man/Grids/#","page":"Manual","title":"Manual","text":"Modules = [Grids]\nPages = [\"grid.jl\"]","category":"page"},{"location":"man/Grids/#Grids.AbstractGrid","page":"Manual","title":"Grids.AbstractGrid","text":"abstract type AbstractGrid{T,N} <: AbstractArray{T,N}\n\nGrids are arrays of points.\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.MappedGrid","page":"Manual","title":"Grids.MappedGrid","text":"A MappedGrid consists of a grid and a map. Each grid point of the mapped grid is the map of the corresponding point of the underlying grid.\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.ProductGrid","page":"Manual","title":"Grids.ProductGrid","text":"A ProductGrid represents the cartesian product of other grids.\n\nstruct ProductGrid{TG,T,N} <: AbstractGrid{T,N}\n\nParameters:\n\nTG is a tuple of (grid) types\nT is the element type of the grid\nN is the dimension of the grid layout\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.ScatteredGrid","page":"Manual","title":"Grids.ScatteredGrid","text":"A grid corresponding to an unstructured collection of points.\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.iscomposite-Tuple{AbstractGrid}","page":"Manual","title":"Grids.iscomposite","text":"iscomposite(grid::AbstractGrid)\n\nDoes the grid consist of multiple grids such as e.g. ProductGrid? Used for pretty printing.\n\n\n\n\n\n","category":"method"},{"location":"man/Grids/#Grids.randomgrid-Tuple{IntervalSets.Domain,Int64}","page":"Manual","title":"Grids.randomgrid","text":"Compute a scattered grid of M points randomly distributed in Ω, using the uniform probability measure on Ω.\n\n\n\n\n\n","category":"method"},{"location":"man/Grids/#Grids.IndexSubGrid","page":"Manual","title":"Grids.IndexSubGrid","text":"struct IndexSubGrid{G,I,T,N} <: AbstractSubGrid{T,N}\n\nAn IndexSubGrid is a subgrid corresponding to a certain range of indices of the underlying grid.\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.MaskedGrid","page":"Manual","title":"Grids.MaskedGrid","text":"struct MaskedGrid{G,M,I,T} <: AbstractSubGrid{T,1}\n\nA MaskedGrid is a subgrid of another grid that is defined by a mask. The mask is true or false for each point in the supergrid. The set of points for which it is true make up the MaskedGrid.\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.hasextension-Tuple{AbstractGrid}","page":"Manual","title":"Grids.hasextension","text":"hasextension(grid::AbstractGrid)\n\nIs it possible to use the resize function. See also resize\n\n\n\n\n\n","category":"method"},{"location":"man/Grids/#Grids.randompoint","page":"Manual","title":"Grids.randompoint","text":"Generate a single random point inside the given domain, with eltype T. Random points are generated inside the given box, until one is inside the domain.\n\n\n\n\n\n","category":"function"},{"location":"man/Grids/#Grids.randompoint-Tuple{DomainSets.ProductDomain}","page":"Manual","title":"Grids.randompoint","text":"Generate a single random point inside the given box, with eltype T.\n\n\n\n\n\n","category":"method"},{"location":"man/Grids/#Grids.resize-Union{Tuple{T}, Tuple{AbstractGrid{T,N} where N,Vararg{Any,N} where N}} where T","page":"Manual","title":"Grids.resize","text":"resize(grid::AbstractGrid, dims...)\n\nCreate a grid of same structure (such as product structure) but with different points in the different dimensions.\n\n\n\n\n\n","category":"method"},{"location":"man/Grids/#","page":"Manual","title":"Manual","text":"Modules = [Grids]\nPages   = [\"intervalgrids.jl\"]","category":"page"},{"location":"man/Grids/#Grids.AbstractEquispacedGrid","page":"Manual","title":"Grids.AbstractEquispacedGrid","text":"abstract type AbstractEquispacedGrid{T} <: AbstractIntervalGrid{T}\n\nAn equispaced grid has equispaced points, and therefore it has a step.\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.AbstractIntervalGrid","page":"Manual","title":"Grids.AbstractIntervalGrid","text":"abstract type AbstractIntervalGrid{T} <: AbstractGrid1d{T}\n\nAn AbstractIntervalGrid is a grid that is defined on an interval, i.e. it is connected.\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.ChebyshevExtremae","page":"Manual","title":"Grids.ChebyshevExtremae","text":"struct ChebyshevExtremae{T} <: AbstractIntervalGrid{T}\n\nA grid with chebyshev extrema on [-1,1].\n\nExample\n\njulia> ChebyshevExtremae(4)\n4-element ChebyshevExtremae{Float64}:\n  1.0\n  0.5000000000000001\n -0.4999999999999998\n -1.0\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.ChebyshevTNodes","page":"Manual","title":"Grids.ChebyshevTNodes","text":"struct ChebyshevNodes{T} <: AbstractIntervalGrid{T}\n\nA grid with chebyshev nodes on [-1,1].\n\nExample\n\njulia> ChebyshevExtremae(4)\n4-element ChebyshevExtremae{Float64}:\n  1.0\n  0.5000000000000001\n -0.4999999999999998\n -1.0\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.EquispacedGrid","page":"Manual","title":"Grids.EquispacedGrid","text":"struct EquispacedGrid{T} <: AbstractEquispacedGrid{T}\n\nAn equispaced grid with n points on an interval [a,b], including the endpoints. It has step (b-a)/(n-1).\n\nExample\n\njulia> EquispacedGrid(4,0,1)\n4-element EquispacedGrid{Float64}:\n 0.0\n 0.3333333333333333\n 0.6666666666666666\n 1.0\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.FourierGrid","page":"Manual","title":"Grids.FourierGrid","text":"struct FourierGrid{T} <: AbstractEquispacedGrid{T}\n\nA Fourier grid is a periodic equispaced grid on the interval [0,1).\n\nexample\n\njulia> FourierGrid(4)\n4-element FourierGrid{Float64}:\n 0.0\n 0.25\n 0.5\n 0.75\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.MidpointEquispacedGrid","page":"Manual","title":"Grids.MidpointEquispacedGrid","text":"struct MidpointEquispacedGrid{T} <: AbstractEquispacedGrid{T}\n\nA MidpointEquispaced grid is an equispaced grid with grid points in the centers of the equispaced subintervals. In other words, this is a DCT-II grid. It has step (b-a)/n.\n\nExample\n\njulia> MidpointEquispacedGrid(4,0,1)\n4-element MidpointEquispacedGrid{Float64}:\n 0.125\n 0.375\n 0.6249999999999999\n 0.875\n\n\n\n\n\n","category":"type"},{"location":"man/Grids/#Grids.PeriodicEquispacedGrid","page":"Manual","title":"Grids.PeriodicEquispacedGrid","text":"struct PeriodicEquispacedGrid{T} <: AbstractEquispacedGrid{T}\n\nA periodic equispaced grid is an equispaced grid that omits the right endpoint. It has step (b-a)/n.\n\nExample\n\njulia> PeriodicEquispacedGrid(4,0,1)\n4-element PeriodicEquispacedGrid{Float64}:\n 0.0\n 0.25\n 0.5\n 0.75\n\n\n\n\n\n","category":"type"}]
}
