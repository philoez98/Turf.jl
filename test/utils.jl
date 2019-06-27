@testset "utils" begin
    @test radiansToLength(1, "radians") == 1
    @test radiansToLength(1, "kilometers") == earthRadius / 1000
    @test radiansToLength(1, "miles") == earthRadius / 1609.344

    @test lengthToRadians(1, "radians") == 1
    @test lengthToRadians(earthRadius / 1000, "kilometers") == 1
    @test lengthToRadians(earthRadius / 1609.344, "miles") == 1

    @test lengthToDegrees(1, "radians") ≈ 57.29 atol=0.01
    @test lengthToDegrees(100, "kilometers") ≈ 0.89 atol=0.01
    @test lengthToDegrees(10, "miles") ≈ 0.14 atol=0.01

    @test convertLength(1000, "meters") == 1
    @test convertLength(1, "kilometers", "miles") ≈ 0.62 atol=0.01
    @test convertLength(1, "nauticalmiles") == 1.852
    @test convertLength(1, "miles", "kilometers") == 1.609344

    @test convertArea(1000) == 0.001
    @test convertArea(1, "kilometers", "miles") == 0.386
    @test convertArea(1, "metres", "centimeters") == 10000
    @test convertArea(100, "meters", "feet") ≈ 1076.39 atol=0.01

    poly = Polygon([[[125, -15],[113, -22],[117, -37],[130, -33],[148, -39],[154, -27],[144, -15],[125, -15]]])
    @test area(poly) ≈ 7766240997209
end
