
# Grids.jl Documentation

For installation instructions, see [Installation](@ref).

For a  full description of the functionality use the manual:
```@contents
Pages = ["man/Grids.md"]
```

## Installation

Grids.jl is not added to the Julia General registry.

### (Recomanded)
For Julia 1.1 or higher, you can add the FrameFun registry and than add Grids.
From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/vincentcp/FrameFunRegistry
pkg> add BasisFunctions
```

### Legacy
In Julia 1.0, the packages can be installed by cloning their git repository. From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/vincentcp/Grids.jl
```

or in a file you could use

```julia
using Pkg
pkg"add https://github.com/vincentcp/Grids.jl"
```
