using Documenter, Turf

makedocs(sitename="Turf.jl", doctest = false, pages=["Home" => ["index.md"],
    "Manual" => ["Getting Started" => "getting-started.md",
    "Examples" => "examples.md",
    hide("points-inside-polygon.md"), hide("simplify.md"),
    hide("nearest-point.md"), hide("grids.md"),
    "Methods" => "methods.md"]], format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"))

deploydocs(
    repo = "github.com/philoez98/Turf.jl.git",
)
