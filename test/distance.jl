include("../src/lib/Distance.jl")
include("../src/geojson/Geometries.jl")


@testset "distance" begin
    p1 = Geometries.Position([-180, -90])
    p2 = Geometries.Position([180, -90])
    @test distance(p1, p1) == 0

    p3 = Geometries.Position([-75.343, 39.984])
    p4 = Geometries.Position([-75.534, 39.123])

    @test distance(p3, p4, "radians") ≈ 0.015 atol=0.001
    @test distance(p3, p4, "degrees") ≈ 0.87 atol=0.01
    @test distance(p3, p4, "miles") ≈ 60.35 atol=0.01

    @test rhumbDistance(p3, p4, "miles") ≈ 60.35 atol=0.01
    @test rhumbDistance(p3, p4, "degrees") ≈ 0.87 atol=0.01
    @test rhumbDistance(p3, p4, "nauticalmiles") ≈ 52.44 atol=0.01

    a =  Geometries.Point([-75.833, 39.284])
    b =  Geometries.Point([-75.6, 39.984])
    c =  Geometries.Point([-75.221, 39.125])
    d =  Geometries.Point([-75.9221, 39.27])
    e =  Geometries.Point([-75.534, 39.123])
    target = Geometries.Point([-75.4, 38.4])

    @test nearestPoint(target, a)

end
