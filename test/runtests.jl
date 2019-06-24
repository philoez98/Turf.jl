using Test, Turf
using GeoInterface: Position, Point, LineString, Polygon, geotype, AbstractGeometry, FeatureCollection

@testset "Turf" begin
    include("angle.jl")
    include("bearing.jl")
    include("utils.jl")
    include("distance.jl")
    include("destination.jl")
    include("booleans.jl")
    include("centering.jl")
end # begin
