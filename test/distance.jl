@testset "distance" begin
    p1 = Position([-180, -90])
    p2 = Position([180, -90])
    @test distance(p1, p1) == 0

    p3 = Position([-75.343, 39.984])
    p4 = Position([-75.534, 39.123])

    @test distance(p3, p4, "radians") ≈ 0.015 atol=0.001
    @test distance(p3, p4, "degrees") ≈ 0.87 atol=0.01
    @test distance(p3, p4, "miles") ≈ 60.35 atol=0.01

    @test rhumb_distance(p3, p4, "miles") ≈ 60.35 atol=0.01
    @test rhumb_distance(p3, p4, "degrees") ≈ 0.87 atol=0.01
    @test rhumb_distance(p3, p4, "nauticalmiles") ≈ 52.44 atol=0.01

    @test rhumb_distance(Point([-75.343, 39.984]), Point([-75.534, 39.123]), "miles") == rhumb_distance(p3, p4, "miles")

    a =  Point([-75.833, 39.284])
    b =  Point([-75.6, 39.984])
    c =  Point([-75.221, 39.125])
    d =  Point([-75.9221, 39.27])
    e =  Point([-75.534, 39.123])
    f = Point([-75.33, 39.44])
    target = Point([-75.4, 39.4])

    @test nearestpoint(target, [a, b, c, d, e, f]).coordinates == [-75.33, 39.44]

    @test pnorm_distance(Point([2, 0]), Point([0, 0]), 2) == 2
    @test pnorm_distance(Point([1, 1]), Point([0, 0]), 1) == 2

    mid = midpoint(Point([0, 10]), Point([0, 0]))

    @test distance(Point([0, 10]), mid) ≈ distance(Point([0, 0]), mid)

    point = Point([0, 0])
    line = LineString([[1, 1],[-1, 1]])

    @test point_to_line_distance(point, line, "miles") ≈ 69.11854715938406 atol=0.05

end
