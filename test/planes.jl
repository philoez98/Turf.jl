@testset "planes" begin
    p = Point([1, 1])
    t1 = Polygon([[[0, 0, 0], [2, 0, 0], [1, 2, 2], [0, 0, 0]]])

    @test planepoint(p, t1) == 1.
end
