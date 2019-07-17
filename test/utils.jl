@testset "utils" begin
    @test radians_to_length(1, "radians") == 1
    @test radians_to_length(1, "kilometers") == earth_radius / 1000
    @test radians_to_length(1, "miles") == earth_radius / 1609.344

    @test length_to_radians(1, "radians") == 1
    @test length_to_radians(earth_radius / 1000, "kilometers") == 1
    @test length_to_radians(earth_radius / 1609.344, "miles") == 1

    @test length_to_degrees(1, "radians") ≈ 57.29 atol=0.01
    @test length_to_degrees(100, "kilometers") ≈ 0.89 atol=0.01
    @test length_to_degrees(10, "miles") ≈ 0.14 atol=0.01

    @test convert_length(1000, "meters") == 1
    @test convert_length(1, "kilometers", "miles") ≈ 0.62 atol=0.01
    @test convert_length(1, "nauticalmiles") == 1.852
    @test convert_length(1, "miles", "kilometers") == 1.609344

    @test convert_area(1000) == 0.001
    @test convert_area(1, "kilometers", "miles") == 0.386
    @test convert_area(1, "metres", "centimeters") == 10000
    @test convert_area(100, "meters", "feet") ≈ 1076.39 atol=0.01

    poly = Polygon([[[125, -15],[113, -22],[117, -37],[130, -33],[148, -39],[154, -27],[144, -15],[125, -15]]])
    @test area(poly) ≈ 7766240997209


    line = LineString([
				[
					0,
					0
				],
				[
					1,1
				],
				[
					2,2
				],
				[
					3,3
				],
				[
					4,4
				],
				[
					4,3
				],
				[
					4,0
				],
				[
					3,
					0
				],
				[
					1,
					0
				],
				[
					0,
					0
				],
				[
					0,
					0
				]
		])

    #@test clean(line).coordinates
end
