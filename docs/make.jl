using Documenter, Homotopy

makedocs(
    format = :html,
    sitename = "Homotopies.jl",
    pages = [
        "Introduction" => "index.md",
        "Interface" => "interface.md",
        "Higher level constructs" => "higherlevelconstructs.md"
        ]
)

deploydocs(
    repo   = "github.com/JuliaHomotopyContinuation/Homotopies.jl.git",
    target = "build",
    julia = "0.6",
    osname = "linux",
    deps   = nothing,
    make   = nothing
)
