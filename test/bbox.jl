@testset "bbox" begin

    line = LineString([[102.0, -10.0], [103.0, 1.0], [104.0, 0.0], [130.0, 4.0]])
    mp = MultiPolygon([[[[102.0, 2.0], [103.0, 2.0], [103.0, 3.0], [102.0, 3.0], [102.0, 2.0]]],
           [[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]], [[100.2, 0.2], [100.8, 0.2], [100.8, 0.8], [100.2, 0.8], [100.2, 0.2]]]])

    p1 = Point([102.0, 0.5])
    poly = Polygon([[[101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0], [101.0, 0.0]]])
    ml = MultiLineString( [[[100.0, 0.0], [101.0, 1.0]],
    [[102.0, 2.0], [103.0, 3.0]]])

    fc = FeatureCollection([Feature(p1), Feature(mp), Feature(poly), Feature(line), Feature(ml)])

    @test bbox(fc) == [100, -10, 130, 4]
    @test bbox(p1) == [102, 0.5, 102, 0.5]
    @test bbox(line) == [102, -10, 130, 4]
    @test bbox(poly) == [100, 0, 101, 1]
    @test bbox(ml) == [100, 0, 103, 3]
    @test bbox(mp) == [100, 0, 103, 3]

    @test bbox(Feature(p1)) == bbox(p1)

    bbpoly = [0, 0, 10, 10]

    coords = bbox_polygon(bbpoly).coordinates
    @test length(coords[1]) === 5
    @test coords[1][1][1] == coords[1][end][1]
    @test coords[1][1][2] == coords[1][1][2]

    @test_throws ErrorException bbox_polygon([1, 2, 3, 4, 5, 6])

end
