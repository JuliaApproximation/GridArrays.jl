# GridArrays.jl

| **Documentation** | **Build Status** | **Coverage** |
|-------------------|------------------|--------------|
| [![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaApproximation.github.io/GridArrays.jl/stable)  [![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaApproximation.github.io/GridArrays.jl/dev) | [![Build Status](https://travis-ci.org/JuliaApproximation/GridArrays.jl.png)](https://travis-ci.org/JuliaApproximation/GridArrays.jl) [![Build status](https://ci.appveyor.com/api/projects/status/gh4ka7m9a7qekqu8?svg=true)](https://ci.appveyor.com/project/JuliaApproximation/GridArrays-jl) | [![Coverage](https://codecov.io/gh/JuliaApproximation/GridArrays.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/GridArrays.jl)  [![Coverage Status](https://coveralls.io/repos/github/JuliaApproximation/GridArrays.jl/badge.svg)](https://coveralls.io/github/JuliaApproximation/GridArrays.jl) |

GridArrays defines a collection of basic grids that act as an array. These
arrays are also associated with a domain as defined by [DomainSets.jl](https://github.com/JuliaApproximation/DomainSets.jl).

The package defines the roots of the classical orthogonal polynomials as arrays,
including `ChebyshevTNodes`, `ChebyshevUNodes`, `LegendreNodes` and others.

```julia
julia> using GridArrays

julia> g1 = EquispacedGrid(5, 0, 1)
5-element EquispacedGrid{Float64}:
 0.0
 0.25
 0.5
 0.75
 1.0

julia> covering(g1)
0.0..1.0

julia> g2 = MidpointEquispacedGrid(5, 0..1)
5-element MidpointEquispacedGrid{Float64}:
 0.1
 0.30000000000000004
 0.5
 0.7000000000000001
 0.9

julia> g3 = ChebyshevNodes(4)
4-element ChebyshevTNodes{Float64}:
 -0.9238795325112867
 -0.38268343236508984
  0.3826834323650897
  0.9238795325112867

julia> covering(g3)
-1.0..1.0 (Chebyshev)
```
