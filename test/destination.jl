using GeoInterface: Position
include("../src/lib/Destination.jl")

@testset "destination" begin

    p1 = Position([-75, 38.10096062273525])
    res1 = destination(p1, 100, 0).coordinates

    p2 = Position([-75, 39])
    res2 = destination(p2, 100, 180).coordinates

    p3 = Position([12, -54])
    res3 = rhumbDestination(p3, -100, 45).coordinates

    p4 = Position([-75, 38.10096062273525])
    res4 = rhumbDestination(p4, 100, 0).coordinates

    @test round.(res1; digits=6) == [-75.0, 39.000281]
    @test round.(res2; digits=5) == [-75, 38.10068]
    @test res3 == [10.90974456038191, -54.63591552764877]
    @test res4 == [-75, 39.00028098645979]
end
