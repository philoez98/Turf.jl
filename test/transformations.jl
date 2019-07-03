@testset "transformations" begin

    p1 = Point([-75.69926351308823,45.43145021122502])
    poly = Polygon([[[0, 29], [3.5, 29], [2.5, 32], [0, 29]]])

    @test transformRotate(geojson=p1, angle=80., pivot=Point([-75.6, 45.3])).coordinates ≈ [-75.433041, 45.391507]
    @test transformRotate(geojson=poly, angle=200., pivot=Point([2.48291015625, 27.819644755099446])).coordinates ≈
        [[[4.331318,25.963562],[1.0811,27.016436],[0.878141,23.89646],[4.331318,25.963562]]] atol = 0.0001

    coll = FeatureCollection([Feature(Point([-75.69926351308823,45.43145021122502])),
            Feature(Polygon([[[0, 29], [3.5, 29], [2.5, 32], [0, 29]]]))])

	@test transformScale(coll, 0.1).features[2].geometry.coordinates ≈ [[[
									1.3495,
									29.675
								],
								[
									1.700667,
									29.675
								],
								[
									1.598959,
									29.975
								],
								[
									1.3495,
									29.675]]]

    @test transformScale(coll, 5., "sw").features[1].geometry.coordinates ≈ [-75.699264,45.43145]

	poly2 = Polygon([[[100, 0], [101, 0], [101, 1], [100, 1], [100, 0]]])
	feats = explode(poly2).features
	res = [Point([100, 0]), Point([101, 0]), Point([101, 1]), Point([100, 1]), Point([100, 0])]


	for i in eachindex(feats)
		@test feats[i].geometry.coordinates == res[i].coordinates
	end

	@test flip(Point([77.34374999999999,43.58039085560784,3000])).coordinates == [43.58039085560784,77.34374999999999,3000]

	p2 = Point([92.46093749999999,54.67383096593114])
	poly2 = Polygon([[[48.1640625,20.632784250388028],[76.640625,20.632784250388028],
            [
              76.640625,
              38.8225909761771
            ],
            [
              48.1640625,
              38.8225909761771
            ],
            [
              48.1640625,
              20.632784250388028]]])

	@test polygonTangents(p2, poly2).features[1].geometry.coordinates ≈ [48.1640625,38.8225909761771]
	@test polygonTangents(p2, poly2).features[2].geometry.coordinates ≈ [76.640625, 20.632784250388028]

	poly3 = Polygon([[[-2.275543, 53.464547 ],
		[-2.275543, 53.489271 ],
		[-2.215118, 53.489271 ],
		[-2.215118, 53.464547 ],
		[-2.275543, 53.464547 ]]])

	@test polygonToLine(poly3).coordinates ≈ [
			[
				-2.275543,
				53.464547
			],
			[
				-2.275543,
				53.489271
			],
			[
				-2.215118,
				53.489271
			],
			[
				-2.215118,
				53.464547
			],
			[
				-2.275543,
				53.464547
			]]

	poly5 = Polygon([[[0, 29], [3.5, 29], [2.5, 32], [0, 29]]])
	poly6 = Polygon([[
        [
          121.46484375,
          -23.88583769986199,
          10
        ],
        [
          125.5078125,
          -27.488781168937983,
          11
        ],
        [
          129.638671875,
          -25.562265014427492,
          12
        ],
        [
          129.2431640625,
          -22.06527806776582,
          13
        ],
        [
          124.93652343749999,
          -20.385825381874263,
          14
        ],
        [
          121.46484375,
          -23.88583769986199,
          10
        ]]])

	@test transformTranslate(poly5, 300, 70).coordinates ≈ [
					[
						[
							2.911836,
							29.922757
						],
						[
							6.411836,
							29.922757
						],
						[
							5.504792,
							32.922757
						],
						[
							2.911836,
							29.922757
						]
					]]

	@test transformTranslate(poly6, 120, 165, 15, false, "miles").coordinates ≈ [
					[
						[
							121.959747,
							-25.563437,
							25
						],
						[
							126.018505,
							-29.166381,
							26
						],
						[
							130.140552,
							-27.239865,
							27
						],
						[
							129.731173,
							-23.742878,
							28
						],
						[
							125.418767,
							-22.063425,
							29
						],
						[
							121.959747,
							-25.563437,
							25
						]
					]]

	poly7 = Polygon([[[
					-20031393.380526427,
					-1860477.2684611906
				],
				[
					19998984.080533516,
					-1906645.233545436
				],
				[
					20014882.98241683,
					-1959233.9090056375
				],
				[
					-20002653.0578912,
					-1937831.5410857894
				],
				[
					-19990117.385252435,
					-1892275.0722278224
				],
				[
					-20031393.380526427,
					-1860477.2684611906
					]]])

	#@test convertTo(poly7, "wgs84").coordinates

	res = convertTo(Point([10, 40]), "mercator")
	@test convertTo(res, "wgs84") ≈ Point([10, 40]) atol=0.5

	pt1 = Point([50, 51])
	pt2 = MultiPoint([[100, 101], [101, 102]])

	fc1 = combine(FeatureCollection([Feature(pt1), Feature(pt2)]))

	@test fc1.features[1].geometry.coordinates ≈ [[50, 51], [100, 101], [101, 102]]

	l1 = LineString([[102.0,-10.0],[130.0,4.0]])
	l2 = LineString([[40.0,-20.0],[150.0,18.0]])

	fc2 = combine(FeatureCollection([Feature(l1), Feature(l2)]))

	@test fc2.features[1].geometry.coordinates == [[[102, -10], [130, 4], [40, -20], [150, 18]]]

	pl1 = Polygon([[
            [20.0, 0.0],
            [101.0, 0.0],
            [101.0, 1.0],
            [100.0, 1.0],
            [100.0, 0.0],
            [20.0, 0.0]]])

	pl2 = MultiPolygon([
		[[
            [30.0, 0.0],
            [102.0, 0.0],
            [103.0, 1.0],
            [30.0, 0.0]
        ]],
        [[
                [20.0, 5.0],
                [92.0, 5.0],
                [93.0, 6.0],
                [20.0, 5.0]
            ],
            [
                [25, 5],
                [30, 5],
                [30, 5.5],
                [25, 5]]]])

	fc3 = combine(FeatureCollection([Feature(pl1), Feature(pl2)]))
	@test fc3.features[1].geometry.coordinates == [[[[20, 0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0], [20, 0],[30.0, 0.0], [102.0, 0.0], [103.0, 1.0], [30.0, 0.0]]]]

	l3 = LineString([[102.0, -10.0],[130.0, 4.0]])
	ml = MultiLineString([[[40.0, -20.0],[150.0, 18.0], [50, -10],[160, 28]]])

	fc4 = combine(FeatureCollection([Feature(l3), Feature(ml)]))
	@test fc4.features[1].geometry.coordinates == [[[102., -10.], [130., 4.], [40., -20.], [150., 18.], [50, -10], [160, 28]]]
end
