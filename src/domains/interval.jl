
"""
    abstract type AbstractIntervalGrid{T} <: AbstractGrid1d{T}

An AbstractIntervalGrid is a grid that is defined on an interval.
"""
abstract type AbstractIntervalGrid{T} <: AbstractGrid1d{T}
end

isperiodic(::AbstractIntervalGrid) = false
covering(grid::AbstractIntervalGrid) = Interval(grid[1], grid[end])
size(grid::AbstractIntervalGrid) = (grid.n,)

"""
    abstract type AbstractEquispacedGrid{T} <: AbstractIntervalGrid{T}

An equispaced grid has equispaced points, and therefore it has a step.
"""
abstract type AbstractEquispacedGrid{T} <: AbstractIntervalGrid{T}
end

range(grid::AbstractEquispacedGrid) = grid.range
==(g1::AbstractEquispacedGrid, g2::AbstractEquispacedGrid) =
    range(g1)==range(g2) && covering(g1) == covering(g2)

size(grid::AbstractEquispacedGrid) = size(range(grid))
step(grid::AbstractEquispacedGrid) = step(range(grid))

@inline function unsafe_grid_getindex(grid::AbstractEquispacedGrid, i::Int)
	@inbounds getindex(range(grid), i)
end

"""
    struct EquispacedGrid{T} <: AbstractEquispacedGrid{T}

An equispaced grid with n points on an interval [a,b], including the endpoints.
It has step (b-a)/(n-1).

# Example
```jldocs
julia> EquispacedGrid(4,0,1)
4-element EquispacedGrid{Float64}:
 0.0
 0.3333333333333333
 0.6666666666666666
 1.0
```
"""
struct EquispacedGrid{T} <: AbstractEquispacedGrid{T}
    # Use StepRangeLen for higher precision
    range   :: LinRange{T}

    EquispacedGrid{T}(n::Int, a, b) where {T} = new(LinRange(T(a),T(b),n))
end

name(g::EquispacedGrid) = "Equispaced grid"


"""
    struct PeriodicEquispacedGrid{T} <: AbstractEquispacedGrid{T}

A periodic equispaced grid is an equispaced grid that omits the right endpoint.
It has step (b-a)/n.

# Example
```jldocs
julia> PeriodicEquispacedGrid(4,0,1)
4-element PeriodicEquispacedGrid{Float64}:
 0.0
 0.25
 0.5
 0.75
```
"""
struct PeriodicEquispacedGrid{T} <: AbstractEquispacedGrid{T}
    range   :: LinRange{T}
    a   ::  T
    b   ::  T

    PeriodicEquispacedGrid{T}(n::Int, a, b) where {T} = new(LinRange(T(a),T(b),n+1)[1:end-1], a, b)
end

name(::PeriodicEquispacedGrid) = "Periodic equispaced grid"
covering(grid::PeriodicEquispacedGrid) = Interval(grid.a, grid.b)
isperiodic(::PeriodicEquispacedGrid) = true

"""
    struct MidpointEquispacedGrid{T} <: AbstractEquispacedGrid{T}

A MidpointEquispaced grid is an equispaced grid with grid points in the centers of the equispaced
subintervals. In other words, this is a DCT-II grid.
It has step `(b-a)/n`.

# Example
```jldocs
julia> MidpointEquispacedGrid(4,0,1)
4-element MidpointEquispacedGrid{Float64}:
 0.125
 0.375
 0.6249999999999999
 0.875
```
"""
struct MidpointEquispacedGrid{T} <: AbstractEquispacedGrid{T}
    range   ::LinRange{T}
    a   ::  T
    b   ::  T

    MidpointEquispacedGrid{T}(n::Int, a, b) where {T} = new(LinRange(T(a),T(b),2n+1)[2:2:end], a, b)
end

name(g::MidpointEquispacedGrid) = "Equispaced midpoints grid"
covering(grid::MidpointEquispacedGrid) = Interval(grid.a, grid.b)
isperiodic(::MidpointEquispacedGrid) = true


"""
    struct FourierGrid{T} <: AbstractEquispacedGrid{T}

A Fourier grid is a periodic equispaced grid on the interval [0,1).

# example
```jldocs
julia> FourierGrid(4)
4-element FourierGrid{Float64}:
 0.0
 0.25
 0.5
 0.75
```
"""
struct FourierGrid{T} <: AbstractEquispacedGrid{T}
    range ::LinRange{T}

    FourierGrid{T}(n::Int) where {T} = new(LinRange(T(0),T(1),n+1)[1:end-1])
end

FourierGrid(n::Int) = FourierGrid{Float64}(n)
similargrid(g::FourierGrid, ::Type{T}, n::Int) where {T} = FourierGrid{T}(n)

FourierGrid(n::Int, d::AbstractInterval) = FourierGrid(n, endpoints(d)...)
FourierGrid(n::Int, a, b) = rescale(FourierGrid{typeof((b-a)/n)}(n), a, b)

name(g::FourierGrid) = "Periodic Fourier grid"
covering(g::FourierGrid{T}) where {T} = UnitInterval{T}()
isperiodic(::FourierGrid) = true



# Grids with flexible support
for GRID in (:PeriodicEquispacedGrid, :MidpointEquispacedGrid, :EquispacedGrid)
    @eval $GRID(n::Int, d::AbstractInterval) =
        $GRID(n, endpoints(d)...)
    @eval similargrid(grid::$GRID, ::Type{T}, n::Int) where {T} =
        $GRID{T}(n, map(T, endpoints(covering(grid)))...)
    @eval rescale(grid::$GRID, a, b) =
        $GRID{promote_type(typeof(a/2),typeof(b/2),eltype(grid))}(length(grid), a, b)
    @eval $GRID(n::Int, a, b) =
        $GRID{promote_type(typeof(a/2),typeof(b/2))}(n, a, b)
    @eval mapped_grid(grid::$GRID, map::AffineMap) =
        $GRID(length(grid), endpoints(map.(covering(grid)))...)
end

# extendable grids
_extension_size(::PeriodicEquispacedGrid, n::Int, factor::Int) = factor*n
_extension_size(::FourierGrid, n::Int, factor::Int) = factor*n
_extension_size(::EquispacedGrid, n::Int, factor::Int) = factor*n-1

for GRID in (:PeriodicEquispacedGrid,:FourierGrid,:EquispacedGrid)
    @eval hasextension(::$GRID) = true
    @eval extend(grid::$GRID, factor::Int) =
        resize(grid, _extension_size(grid, length(grid), factor))
end

# function mapped_grid(grid::FourierGrid{T}, map::AffineMap) where T
#     s = map*covering(grid)
#     s≈UnitInterval{T}() ?
#         grid : PeriodicEquispacedGrid{T}(length(grid), endpoints(s)...)
# end

function rescale(g::FourierGrid, a, b)
	m = mapto(covering(g), a..b)
	mapped_grid(g, m)
end

mapped_grid(g::FourierGrid, map::AffineMap) = MappedGrid(g, map)
