using Pkg
Pkg.develop(PackageSpec(path=splitdir(@__DIR__)[1]))
pkg"instantiate"
using Documenter, Grids

const format = Documenter.HTML(
        prettyurls = ("deploy" in ARGS),
    )

makedocs(sitename="Grids.jl",
    modules = [Grids],
    format = format,
    pages = [
        "Home" => "index.md",
        "Manual" => "man/Grids.md"
        ],
    doctest=true
)

if "deploy" in ARGS && Sys.ARCH === :x86_64 && Sys.KERNEL === :Linux
    deploydocs(
        repo = "github.com/vincentcp/Grids.jl.git",
    )
end
