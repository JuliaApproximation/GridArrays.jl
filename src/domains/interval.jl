
"""
    abstract type AbstractIntervalGrid{T} <: AbstractGrid1d{T}

An AbstractIntervalGrid is a grid that is defined on an interval.
"""
abstract type AbstractIntervalGrid{T} <: AbstractGrid1d{T} end

isperiodic(::AbstractIntervalGrid) = false


"""
    abstract type AbstractEquispacedGrid{T} <: AbstractIntervalGrid{T}

An equispaced grid has equispaced points, and therefore it has a step.
"""
abstract type AbstractEquispacedGrid{T} <: AbstractIntervalGrid{T} end

==(g1::AbstractEquispacedGrid, g2::AbstractEquispacedGrid) =
    covering(g1) == covering(g2) && size(g1)==size(g2) && (g1[1]==g2[1])


# Mainly for technical reasons we introduce two intermediate types
abstract type AbstractEquispacedRangeGrid{T} <: AbstractEquispacedGrid{T} end
abstract type AbstractUnitEquispacedGrid{T} <: AbstractEquispacedGrid{T} end

# the range grids are defined in terms of a range
range(g::AbstractEquispacedRangeGrid) = g.range
size(g::AbstractEquispacedRangeGrid) = size(range(g))
step(g::AbstractEquispacedRangeGrid) = step(range(g))

@inline function unsafe_grid_getindex(grid::AbstractEquispacedRangeGrid, i::Int)
	@inbounds getindex(range(grid), i)
end

# the unit equispaced grids are defined on the unit interval [0,1]
covering(g::AbstractUnitEquispacedGrid{T}) where {T} = UnitInterval{T}()
size(g::AbstractUnitEquispacedGrid) = (g.n,)

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
struct EquispacedGrid{T} <: AbstractEquispacedRangeGrid{T}
    # Use StepRangeLen for higher precision
    range   :: LinRange{T}

    EquispacedGrid{T}(n::Int, a, b) where {T} = new(LinRange(T(a),T(b),n))
end

covering(g::EquispacedGrid) = g[1]..g[end]

name(g::EquispacedGrid) = "Equispaced grid"


"An equispaced grid on the unit interval `[0,1]`."
struct UnitEquispacedGrid{T} <: AbstractUnitEquispacedGrid{T}
	n	::	Int
end

UnitEquispacedGrid(n::Int) = UnitEquispacedGrid{Float64}(n)
similargrid(g::UnitEquispacedGrid, ::Type{T}, n::Int) where {T} =
	UnitEquispacedGrid{T}(n)

unsafe_grid_getindex(grid::UnitEquispacedGrid{T}, i::Int) where {T} =
	convert(T, (i-1))/(grid.n-1)

"""
    struct PeriodicEquispacedGrid{T} <: AbstractEquispacedGrid{T}

A periodic equispaced grid is an equispaced grid that omits the right endpoint.
It has step size (b-a)/n.

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
struct PeriodicEquispacedGrid{T} <: AbstractEquispacedRangeGrid{T}
    range   :: LinRange{T}
    a   ::  T
    b   ::  T

    PeriodicEquispacedGrid{T}(n::Int, a, b) where {T} =
		new(LinRange(T(a),T(b),n+1)[1:end-1], a, b)
end

range(g::PeriodicEquispacedGrid) = g.range

name(::PeriodicEquispacedGrid) = "Periodic equispaced grid"
covering(grid::PeriodicEquispacedGrid) = grid.a..grid.b
isperiodic(::PeriodicEquispacedGrid) = true


"A periodic equispaced grid on the unit interval `[0,1]`."
struct UnitPeriodicEquispacedGrid{T} <: AbstractUnitEquispacedGrid{T}
	n	::	Int
end

UnitPeriodicEquispacedGrid(n::Int) = UnitPeriodicEquispacedGrid{Float64}(n)
similargrid(g::UnitPeriodicEquispacedGrid, ::Type{T}, n::Int) where {T} =
	UnitPeriodicEquispacedGrid{T}(n)
isperiodic(::UnitPeriodicEquispacedGrid) = true

unsafe_grid_getindex(grid::UnitPeriodicEquispacedGrid{T}, i::Int) where {T} =
	convert(T, (i-1))/grid.n

const FourierGrid = UnitPeriodicEquispacedGrid

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
struct MidpointEquispacedGrid{T} <: AbstractEquispacedRangeGrid{T}
    range   ::LinRange{T}
    a   ::  T
    b   ::  T

    MidpointEquispacedGrid{T}(n::Int, a, b) where {T} =
		new(LinRange(T(a),T(b),2n+1)[2:2:end], a, b)
end

range(g::MidpointEquispacedGrid) = g.range

name(g::MidpointEquispacedGrid) = "Equispaced midpoints grid"
covering(grid::MidpointEquispacedGrid) = grid.a..grid.b
isperiodic(::MidpointEquispacedGrid) = true


"A midpoint equispaced grid on the unit interval `[0,1]`."
struct UnitMidpointEquispacedGrid{T} <: AbstractUnitEquispacedGrid{T}
	n	::	Int
end

UnitMidpointEquispacedGrid(n::Int) = UnitMidpointEquispacedGrid{Float64}(n)
similargrid(g::UnitMidpointEquispacedGrid, ::Type{T}, n::Int) where {T} = UnitMidpointEquispacedGrid{T}(n)

unsafe_grid_getindex(grid::UnitMidpointEquispacedGrid{T}, i::Int) where {T} =
	convert(T, 2(i-1))/(2grid.n-1)


# Grids with flexible support
for GRID in (:PeriodicEquispacedGrid, :MidpointEquispacedGrid, :EquispacedGrid)
	@eval $GRID(n::Int, a, b) = $GRID(n, promote(a, b)...)
	@eval $GRID(n::Int, a::T, b::T) where {T} = $GRID{float(T)}(n, a, b)
	@eval $GRID(n::Int, d::AbstractInterval) =
		$GRID(n, infimum(d), supremum(d))
	@eval $GRID{T}(n::Int, d::AbstractInterval) where {T} =
		$GRID{T}(n, infimum(d), supremum(d))
    @eval similargrid(grid::$GRID, ::Type{T}, n::Int) where {T} =
        $GRID{T}(n, covering(grid))
    @eval rescale(grid::$GRID{T}, a::S, b::U) where {U,S,T} =
        $GRID{promote_type(float(S),float(U),T)}(length(grid), a, b)
    @eval map_grid(grid::$GRID, map) =
        $GRID(length(grid), map_domain(map, covering(grid)))
end

map_grid(g::AbstractEquispacedRangeGrid, map) =
	rescale(g, endpoints(map_domain(map, covering(g)))...)

# extensible grids
# TODO: deprecate extend (in favour of resize)
extend(grid::AbstractEquispacedRangeGrid, factor::Int) =
	resize(grid, extension_size(grid, length(grid), factor))

hasextension(::Union{PeriodicEquispacedGrid,FourierGrid,EquispacedGrid}) = true

extension_size(::PeriodicEquispacedGrid, n::Int, factor::Int) = factor*n
extension_size(::FourierGrid, n::Int, factor::Int) = factor*n
extension_size(::EquispacedGrid, n::Int, factor::Int) = factor*n-1
