
"""
    abstract type AbstractIntervalGrid{T} <: AbstractGrid1d{T}

An AbstractIntervalGrid is a grid that is defined on an interval, i.e. it is connected.
"""
abstract type AbstractIntervalGrid{T} <: AbstractGrid1d{T}
end

isperiodic(::AbstractIntervalGrid) = false
support(grid::AbstractIntervalGrid) = Interval(grid[1], grid[end])
size(grid::AbstractIntervalGrid) = (grid.n,)

string(g::AbstractIntervalGrid) = name(g) * " of length $(length(g)) on [$(support(g)))], ELT = $(eltype(g))"


instantiate(::Type{T}, n::Int, ::Type{ELT})  where {T<:AbstractIntervalGrid,ELT} = T(n,ELT(0),ELT(1))

"""
    abstract type AbstractEquispacedGrid{T} <: AbstractIntervalGrid{T}

An equispaced grid has equispaced points, and therefore it has a step.
"""
abstract type AbstractEquispacedGrid{T} <: AbstractIntervalGrid{T}
end

range(grid::AbstractEquispacedGrid) = grid.range
==(g1::AbstractEquispacedGrid, g2::AbstractEquispacedGrid) =
    range(g1)==range(g2) && support(g1) == support(g2)

size(grid::AbstractEquispacedGrid) = (length(range(grid)),)
step(grid::AbstractEquispacedGrid) = step(range(grid))
unsafe_getindex(grid::AbstractEquispacedGrid, i::Int) = unsafe_getindex(range(grid), i)

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
support(grid::PeriodicEquispacedGrid) = Interval(grid.a, grid.b)
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
support(grid::MidpointEquispacedGrid) = Interval(grid.a, grid.b)
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

name(g::FourierGrid) = "Periodic Fourier grid"
support(g::FourierGrid{T}) where {T} = UnitInterval{T}()
isperiodic(::FourierGrid) = true

"""
    struct ChebyshevNodes{T} <: AbstractIntervalGrid{T}

A grid with chebyshev nodes on [-1,1].

# Example
```jldocs
julia> ChebyshevExtremae(4)
4-element ChebyshevExtremae{Float64}:
  1.0
  0.5000000000000001
 -0.4999999999999998
 -1.0
```
"""
struct ChebyshevNodes{T} <: AbstractIntervalGrid{T}
    n   ::  Int
end

const ChebyshevGrid = ChebyshevNodes
const ChebyshevPoints = ChebyshevNodes

support(::ChebyshevNodes{T}) where T = ChebyshevInterval{T}()

# The minus sign is added to avoid having to flip the inputs to the dct. More elegant fix required.
unsafe_getindex(g::ChebyshevNodes{T}, i::Int) where {T} = T(-1)*cos((i-T(1)/2) * T(pi) / (g.n) )
name(g::ChebyshevNodes) = "Chebyshev nodes"


"""
    struct ChebyshevExtremae{T} <: AbstractIntervalGrid{T}

A grid with chebyshev extrema on [-1,1].

# Example
```jldocs
julia> ChebyshevExtremae(4)
4-element ChebyshevExtremae{Float64}:
  1.0
  0.5000000000000001
 -0.4999999999999998
 -1.0
```
"""
struct ChebyshevExtremae{T} <: AbstractIntervalGrid{T}
    n   ::  Int
end

const ChebyshevPointsOfTheSecondKind = ChebyshevExtremae

name(g::ChebyshevExtremae) = "Chebyshev extremae"
minimum(g::ChebyshevExtremae{T}) where {T} = -one(T)
maximum(g::ChebyshevExtremae{T}) where {T} = one(T)
support(::ChebyshevExtremae{T}) where T = ChebyshevInterval{T}()

# TODO: flip the values so that they are sorted
unsafe_getindex(g::ChebyshevExtremae{T}, i::Int) where {T} = i == 0 ? T(0) : cos((i-1)*T(pi) / (g.n-1) )

# Grids with flexible support
for GRID in (:PeriodicEquispacedGrid, :MidpointEquispacedGrid, :EquispacedGrid)
    @eval $GRID(n::Int, d::AbstractInterval{T}) where {T} =
        $GRID(n, endpoints(d)...)
    @eval similargrid(grid::$GRID, ::Type{T}, n::Int) where {T} =
        $GRID{T}(n, map(T, endpoints(support(grid)))...)
    @eval rescale(grid::$GRID, a, b) =
        $GRID{promote_type(typeof(a/2),typeof(b/2),eltype(grid))}(length(grid), a, b)
    @eval $GRID(n::Int, a, b) =
        $GRID{promote_type(typeof(a/2),typeof(b/2))}(n, a, b)
    @eval mapped_grid(grid::$GRID, map::AffineMap) =
        $GRID(length(grid), endpoints(map*support(grid))...)
end

# Grids with fixed support
for GRID in (:FourierGrid, :ChebyshevNodes, :ChebyshevExtremae)
    @eval similargrid(g::$GRID, ::Type{T}, n::Int) where {T} = $GRID{T}(n)
    @eval $GRID(n::Int) = $GRID{Float64}(n)
    @eval $GRID(n::Int, a, b) = rescale($GRID{typeof((b-a)/n)}(n), a, b)
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

function mapped_grid(grid::FourierGrid{T}, map::AffineMap) where T
    s = map*support(grid)
    s≈UnitInterval{T}() ?
        grid : PeriodicEquispacedGrid{T}(length(grid), endpoints(s)...)
end
