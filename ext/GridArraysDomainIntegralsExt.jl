module GridArraysDomainIntegralsExt

using GridArrays, DomainIntegrals

import DomainIntegrals:
    measure,
    laguerre_α,
    jacobi_α,
    jacobi_β

measure(x::GridArrays.LaguerreNodes) = LaguerreWeight(x.α)
measure(w::GridArrays.LaguerreWeights) = LaguerreWeight(w.α)

laguerre_α(x::GridArrays.LaguerreNodes) = x.α
laguerre_α(w::GridArrays.LaguerreWeights) = w.α

measure(x::GridArrays.JacobiNodes) = JacobiWeight(x.α, x.β)
measure(w::GridArrays.JacobiWeights) = JacobiWeight(w.α, w.β)

jacobi_α(x::GridArrays.JacobiNodes) = x.α
jacobi_α(w::GridArrays.JacobiWeights) = w.α
jacobi_β(x::GridArrays.JacobiNodes) = x.β
jacobi_β(w::GridArrays.JacobiWeights) = w.β

measure(x::GridArrays.LegendreNodes{T}) where T = LegendreWeight{T}()
measure(w::GridArrays.LegendreWeights{T}) where T = LegendreWeight{T}()

jacobi_α(x::GridArrays.LegendreNodes{T}) where T = zero(T)
jacobi_α(w::GridArrays.LegendreWeights{T}) where T = zero(T)
jacobi_β(x::GridArrays.LegendreNodes{T}) where T = zero(T)
jacobi_β(w::GridArrays.LegendreWeights{T}) where T = zero(T)

measure(x::GridArrays.ChebyshevTNodes{T}) where T = ChebyshevTWeight{T}()
measure(w::GridArrays.ChebyshevTWeights{T}) where T = ChebyshevTWeight{T}()

jacobi_α(x::GridArrays.ChebyshevTNodes{T}) where T = -one(T)/2
jacobi_α(w::GridArrays.ChebyshevTWeights{T}) where T = -one(T)/2
jacobi_β(x::GridArrays.ChebyshevTNodes{T}) where T = -one(T)/2
jacobi_β(w::GridArrays.ChebyshevTWeights{T}) where T = -one(T)/2

measure(x::GridArrays.ChebyshevUNodes{T}) where T = ChebyshevUWeight{T}()
measure(w::GridArrays.ChebyshevUWeights{T}) where T = ChebyshevUWeight{T}()

jacobi_α(x::GridArrays.ChebyshevUNodes{T}) where T = -one(T)/2
jacobi_α(w::GridArrays.ChebyshevUWeights{T}) where T = -one(T)/2
jacobi_β(x::GridArrays.ChebyshevUNodes{T}) where T = -one(T)/2
jacobi_β(w::GridArrays.ChebyshevUWeights{T}) where T = -one(T)/2

measure(x::GridArrays.HermiteNodes{T}) where T = HermiteWeight{T}()
measure(w::GridArrays.HermiteWeights{T}) where T = HermiteWeight{T}()

end
