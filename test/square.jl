@testset "square" begin
    b1 = [0, 0, 5, 10]
    b2 = [0, 0, 10, 5]

    @test square(b1) == [-2.5, 0, 7.5, 10]
    @test square(b2) == [0, -2.5, 10, 7.5]
end
