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

    @test_throws ErrorException convert_area(-5, "meters", "feet")

    poly = Polygon([[[125, -15],[113, -22],[117, -37],[130, -33],[148, -39],[154, -27],[144, -15],[125, -15]]])
    @test area(poly) ≈ 7766240997209


    line = LineString([
	[
		0,
		0
	],
	[
		0,
		2
	],
	[
		0,
		5
	],
	[
		0,
		8
	],
	[
		0,
		10
	],
	[
		0,
		10
		]])

	poly = Polygon([
			[
				[
					0,
					0
				],
				[
					0,
					5
				],
				[
					0,
					10
				],
				[
					10,
					10
				],
				[
					10,
					0
				],
				[
					5,
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
			]
		])

	ml = MultiLineString([
			[
				[
					0,
					0
				],
				[
					0,
					2
				],
				[
					0,
					5
				],
				[
					0,
					8
				],
				[
					0,
					10
				],
				[
					0,
					10
				]
			],
			[
				[
					1,
					1
				],
				[
					2,
					2
				],
				[
					3,
					3
				],
				[
					4,
					4
				],
				[
					5,
					5
				],
				[
					6,
					6
				]
			]
		])

	mp = MultiPolygon([
			[
				[
					[
						0,
						0
					],
					[
						5,
						0
					],
					[
						5,
						0
					],
					[
						10,
						0
					],
					[
						10,
						10
					],
					[
						0,
						10
					],
					[
						0,
						5
					],
					[
						0,
						0
					]
				],
				[
					[
						1,
						5
					],
					[
						1,
						7
					],
					[
						1,
						8.5
					],
					[
						4.5,
						8.5
					],
					[
						4.5,
						7
					],
					[
						4.5,
						5
					],
					[
						1,
						5
					]
				]
			],
			[
				[
					[
						11,
						11
					],
					[
						11.5,
						11.5
					],
					[
						12,
						12
					],
					[
						12,
						11
					],
					[
						11.5,
						11
					],
					[
						11,
						11
					],
					[
						11,
						11
					]
				]
			]
		])

	mpo = MultiPoint([
			[
				14.765625,
				26.194876675795218
			],
			[
				8.61328125,
				23.483400654325642
			],
			[
				17.75390625,
				24.926294766395593
			]
		])
    @test clean(line).coordinates == [[0, 0], [0, 10]]
	@test clean(poly).coordinates == [
			[
				[
					0,
					0
				],
				[
					0,
					10
				],
				[
					10,
					10
				],
				[
					10,
					0
				],
				[
					0,
					0
				]
			]
		]

	@test clean(ml).coordinates == [
			[
				[
					0,
					0
				],
				[
					0,
					10
				]
			],
			[
				[
					1,
					1
				],
				[
					6,
					6
				]
			]
		]
	@test clean(mp).coordinates == [
			[
				[
					[
						0,
						0
					],
					[
						10,
						0
					],
					[
						10,
						10
					],
					[
						0,
						10
					],
					[
						0,
						0
					]
				],
				[
					[
						1,
						5
					],
					[
						1,
						8.5
					],
					[
						4.5,
						8.5
					],
					[
						4.5,
						5
					],
					[
						1,
						5
					]
			],
			[
					[
						11,
						11
					],
					[
						12,
						12
					],
					[
						12,
						11
					],
					[
						11,
						11
					]
			]
		]
	]

	@test clean(mpo).coordinates == mpo.coordinates
	@test clean(Feature(mpo)).coordinates == clean(mpo).coordinates
end
