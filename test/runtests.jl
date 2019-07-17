using Test, Turf

@testset "Turf" begin
    include("angle.jl")
    include("bearing.jl")
    include("utils.jl")
    include("distance.jl")
    include("destination.jl")
    include("booleans.jl")
    include("centering.jl")
    include("transformations.jl")
    include("ellipse.jl")
    include("circle.jl")
    include("planes.jl")
    include("grids.jl")
    include("lines.jl")
    #include("bezier.jl")
    include("square.jl")
end # begin
