using Pkg
Pkg.develop(PackageSpec(path=splitdir(@__DIR__)[1]))
pkg"instantiate"
using Documenter, Grids

const render_pdf = "pdf" in ARGS
let r = r"buildroot=(.+)", i = findfirst(x -> occursin(r, x), ARGS)
    global const buildroot = i === nothing ? (@__DIR__) : first(match(r, ARGS[i]).captures)
end

const format = if render_pdf
    LaTeX(
        platform = "texplatform=docker" in ARGS ? "docker" : "native"
    )
else
    Documenter.HTML(
        prettyurls = ("deploy" in ARGS),
    )
end

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
