using Documenter, Turf

makedocs(sitename="Turf.jl", pages=["Home" => ["index.md"],
    "Manual" => ["Getting Started" => "getting-started.md",
    "Examples" => "examples.md",
    "Methods" => "methods.md"]], format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"))

deploydocs(
    repo = "github.com/philoez98/Turf.jl.git",
)
