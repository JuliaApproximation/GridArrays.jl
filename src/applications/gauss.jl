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
struct ChebyshevTNodes{T} <: AbstractIntervalGrid{T}
    n   ::  Int
end

const ChebyshevNodes = ChebyshevTNodes

size(g::ChebyshevTNodes) = (g.n,)

unsafe_grid_getindex(g::ChebyshevTNodes{T}, i::Int) where {T} = T(-1)*cos((i-T(1)/2) * T(pi) / (g.n) )


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

size(g::ChebyshevExtremae) = (g.n,)

# TODO: flip the values so that they are sorted
unsafe_grid_getindex(g::ChebyshevExtremae{T}, i::Int) where {T} = i == 0 ? T(0) : cos((i-1)*T(pi) / (g.n-1) )


struct ChebyshevUNodes{T} <: AbstractIntervalGrid{T}
    n   :: Int
end

size(g::ChebyshevUNodes) = (g.n,)

unsafe_grid_getindex(nodes::ChebyshevUNodes{T}, i::Int) where T = cos((nodes.n + 1 - i) * convert(T,π) / (nodes.n + 1))


struct LegendreNodes{T} <: AbstractIntervalGrid{T}
    nodes   :: Vector{T}
end

size(grid::LegendreNodes) = size(grid.nodes)
unsafe_grid_getindex(grid::LegendreNodes, i::Int) = @inbounds getindex(grid.nodes, i)


struct LaguerreNodes{T} <: AbstractIntervalGrid{T}
    α       :: T
    nodes   :: Vector{T}
end

size(grid::LaguerreNodes) = size(grid.nodes)
unsafe_grid_getindex(grid::LaguerreNodes, i::Int) = @inbounds getindex(grid.nodes, i)


struct HermiteNodes{T} <: AbstractIntervalGrid{T}
    nodes   :: Vector{T}
end

size(grid::HermiteNodes) = size(grid.nodes)
unsafe_grid_getindex(grid::HermiteNodes, i::Int) = @inbounds getindex(grid.nodes, i)


struct JacobiNodes{T} <: AbstractIntervalGrid{T}
    α   ::  T
    β   ::  T
    nodes   ::Vector
end

size(grid::JacobiNodes) = size(grid.nodes)
unsafe_grid_getindex(grid::JacobiNodes, i::Int) = @inbounds getindex(grid.nodes, i)



"A vector that aliases another vector but has its own type."
abstract type VectorAlias{T} <: AbstractVector{T} end
@propagate_inbounds getindex(vector::VectorAlias, i::Int) = getindex(vector.vector, i)
size(vector::VectorAlias) = size(vector.vector)

"A function vector is a lazy vector with entries defined by a function."
abstract type FunctionVector{T} <: AbstractVector{T} end
size(vector::FunctionVector) = (vector.n,)

@propagate_inbounds function getindex(vector::FunctionVector, i::Int)
    @boundscheck (1 <= i <= length(vector)) || throw(BoundsError())
    funvector_getindex(vector, i)
end


struct ChebyshevTWeights{T,F} <: VectorAlias{T}
    vector  ::  F
end

# The Chebyshev weights are constant
function ChebyshevTWeights{T}(n) where {T}
    vector = Fill(convert(T,π) / n, n)
    ChebyshevTWeights{T,typeof(vector)}(vector)
end

struct ChebyshevUWeights{T} <: FunctionVector{T}
    n   :: Int
end

funvector_getindex(weights::ChebyshevUWeights{T}, i::Int) where {T} =
    convert(T,π)/(weights.n + 1) * sin(convert(T,weights.n + 1 -i) / (weights.n + 1) * convert(T,π))^2

gausschebyshev(::Type{T}, n::Int) where T =
    ChebyshevTNodes{T}(n), ChebyshevTWeights{T}(n)
gausschebyshev(n::Int) = gausschebyshev(Float64, n)


gausschebyshevu(::Type{T}, n::Int) where T =
    ChebyshevUNodes{T}(n), ChebyshevUWeights{T}(n)
gausschebyshevu(n::Int) = gausschebyshevu(Float64,n)

struct LegendreWeights{T} <: VectorAlias{T}
    vector  ::  Vector{T}
end

function gausslegendre(::Type{Float64}, n::Int)
    x,w = FastGaussQuadrature.gausslegendre(n)
    LegendreNodes(x), LegendreWeights(w)
end
function gausslegendre(::Type{T}, n::Int) where T
    x,w = GaussQuadrature.legendre(T, n)
    LegendreNodes(x), LegendreWeights(w)
end
gausslegendre(n::Int) = gausslegendre(Float64, n)

struct LaguerreWeights{T} <: VectorAlias{T}
    α       ::  T
    vector  ::  Vector{T}
end

function gausslaguerre(::Type{Float64}, n::Int, α::Float64)
    x,w = FastGaussQuadrature.gausslaguerre(n, α)
    LaguerreNodes(α, x), LaguerreWeights(α, w)
end
function gausslaguerre(::Type{T}, n::Int, α::T) where T
    x,w = GaussQuadrature.laguerre(n, α)
    LaguerreNodes(α, x), LaguerreWeights(α, w)
end
gausslaguerre(n::Int, α::T) where T = gausslaguerre(T, n, α)

struct HermiteWeights{T} <: VectorAlias{T}
    vector  ::  Vector{T}
end
function gausshermite(::Type{Float64}, n::Int)
    x,w = FastGaussQuadrature.gausshermite(n)
    HermiteNodes(x), HermiteWeights(w)
end
function gausshermite(::Type{T}, n::Int) where T
    x,w = GaussQuadrature.hermite(T, n)
    HermiteNodes(x), HermiteWeights(w)
end
gausshermite(n::Int) = gausshermite(Float64, n)

struct JacobiWeights{T} <: VectorAlias{T}
    α       ::  T
    β       ::  T
    vector  ::  Vector{T}
end
function gaussjacobi(::Type{Float64}, n::Int, α::Float64, β::Float64)
    x,w = FastGaussQuadrature.gaussjacobi(n, α, β)
    JacobiNodes(α, β, x), JacobiWeights(α, β, w)
end
function gaussjacobi(::Type{T}, n::Int, α::T, β::T) where T
    x,w = GaussQuadrature.jacobi(n, α, β)
    JacobiNodes(α, β, x), JacobiWeights(α, β, w)
end
gaussjacobi(n::Int, α::T, β::T) where T = gaussjacobi(T, n, α, β)



covering(::ChebyshevTNodes{T}) where T = ChebyshevInterval{T}()
name(g::ChebyshevTNodes) = "ChebyshevT nodes"

name(g::ChebyshevExtremae) = "Chebyshev extremae"
covering(::ChebyshevExtremae{T}) where T = ChebyshevInterval{T}()

name(g::ChebyshevUNodes) = "ChebyshevU nodes"
covering(::ChebyshevUNodes{T}) where T = ChebyshevInterval{T}()

name(g::LegendreNodes) = "Legendre nodes"
covering(::LegendreNodes{T}) where T = ChebyshevInterval{T}()
LegendreNodes{T}(n::Int) where T = gausslegendre(T, n)[1]

name(g::LaguerreNodes) = "Laguerre nodes α=$(g.α)"
covering(::LaguerreNodes{T}) where T = HalfLine{T}()
similargrid(grid::LaguerreNodes, T, n::Int) = LaguerreNodes{T}(n, T(grid.α))
LaguerreNodes(n::Int, α::T) where T = gausslaguerre(T, n, α)[1]
LaguerreNodes{T}(n::Int, α::T) where T = gausslaguerre(T, n, α)[1]

name(g::HermiteNodes) = "Hermite nodes"
covering(::HermiteNodes{T}) where T = DomainSets.FullSpace{T}()
HermiteNodes{T}(n::Int) where T = gausshermite(T, n)[1]

name(g::JacobiNodes) = "Jacobi nodes α=$(g.α), β=$(g.β)"
covering(::JacobiNodes{T}) where T = ChebyshevInterval{T}()
similargrid(grid::JacobiNodes, T, n::Int) = JacobiNodes{T}(n, T(grid.α), T(grid.β))
JacobiNodes{T}(n::Int, α::T, β::T) where T = gaussjacobi(T, n, α, β)[1]
JacobiNodes(n::Int, α::T, β::T) where T = gaussjacobi(T, n, α, β)[1]


# Grids with fixed support and one variable
for GRID in (:ChebyshevNodes, :ChebyshevExtremae, :ChebyshevUNodes, :LegendreNodes, :HermiteNodes)
    @eval similargrid(g::$GRID, ::Type{T}, n::Int) where {T} = $GRID{T}(n)
    @eval $GRID(n::Int) = $GRID{Float64}(n)
    # @eval $GRID(n::Int, d::AbstractInterval) =
        # $GRID(n, endpoints(d)...)
    # @eval $GRID(n::Int, a, b) = rescale($GRID{typeof((b-a)/n)}(n), a, b)
end
