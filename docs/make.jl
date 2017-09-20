using Documenter, Homotopy

makedocs(
    format = :html,
    sitename = "Homotopy.jl",
    pages = [
        "Introduction" => "index.md",
        "Interface" => "interface.md",
        ]
)

deploydocs(
    repo   = "github.com/JuliaHomotopyContinuation/Homotopy.jl.git",
    target = "build",
    julia = "0.6",
    osname = "linux",
    deps   = nothing,
    make   = nothing
)
