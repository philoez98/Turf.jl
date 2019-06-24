using Test, Turf
using GeoInterface

@testset "Turf" begin
    include("angle.jl")
    include("bearing.jl")
    include("utils.jl")
    include("distance.jl")
    include("destination.jl")
    include("booleans.jl")
    include("centering.jl")
    include("transformations.jl")
    include("circle.jl")
end # begin
