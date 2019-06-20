include("../src/Constants.jl")
include("../src/Utils.jl")

@testset "utils" begin
    @test radiansToLength(1, "radians") == 1
    @test radiansToLength(1, "kilometers") == Constants.earthRadius / 1000
    @test radiansToLength(1, "miles") == Constants.earthRadius / 1609.344

    @test lengthToRadians(1, "radians") == 1
    @test lengthToRadians(Constants.earthRadius / 1000, "kilometers") == 1
    @test lengthToRadians(Constants.earthRadius / 1609.344, "miles") == 1

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
end
