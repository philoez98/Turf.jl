@testset "distance" begin
    p1 = Position([-180, -90])
    p2 = Position([180, -90])
    @test distance(p1, p1) == 0

    p3 = Position([-75.343, 39.984])
    p4 = Position([-75.534, 39.123])

    @test distance(p3, p4, "radians") ≈ 0.015 atol=0.001
    @test distance(p3, p4, "degrees") ≈ 0.87 atol=0.01
    @test distance(p3, p4, "miles") ≈ 60.35 atol=0.01

    @test rhumbDistance(p3, p4, "miles") ≈ 60.35 atol=0.01
    @test rhumbDistance(p3, p4, "degrees") ≈ 0.87 atol=0.01
    @test rhumbDistance(p3, p4, "nauticalmiles") ≈ 52.44 atol=0.01

    a =  Point([-75.833, 39.284])
    b =  Point([-75.6, 39.984])
    c =  Point([-75.221, 39.125])
    d =  Point([-75.9221, 39.27])
    e =  Point([-75.534, 39.123])
    f = Point([-75.33, 39.44])
    target = Point([-75.4, 39.4])

    @test nearestPoint(target, [a, b, c, d, e, f]).coordinates == [-75.33, 39.44]

    @test pNormDistance(Point([2, 0]), Point([0, 0]), 2) == 2
    @test pNormDistance(Point([1, 1]), Point([0, 0]), 1) == 2

end
