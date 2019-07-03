@testset "angle" begin
    p1 = Position([5., 5.])
    p2 = Position([5., 6.])
    p3 = Position([3., 4.])

    @test round(angle_adjacent(p1, p2, p3, false, false)) == 45.0
    @test round(angle_adjacent(p3, p2, p1, true, false)) == 360 - 45.0
end
