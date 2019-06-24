@testset "transformations" begin

    p1 = Point([-75.69926351308823,45.43145021122502])
    poly = Polygon([[[0, 29], [3.5, 29], [2.5, 32], [0, 29]]])

    @test transformRotate(geojson=p1, angle=80., pivot=Point([-75.6, 45.3])).coordinates ≈ [-75.433041, 45.391507]
    @test transformRotate(geojson=poly, angle=200., pivot=Point([2.48291015625, 27.819644755099446])).coordinates ≈
        [[[4.331318,25.963562],[1.0811,27.016436],[0.878141,23.89646],[4.331318,25.963562]]] atol = 0.0001

    coll = FeatureCollection([Feature(Point([-75.69926351308823,45.43145021122502])),
            Feature(Polygon([[[0, 29], [3.5, 29], [2.5, 32], [0, 29]]]))])

    @test transformScale(coll, 5., "sw").features[1].geometry.coordinates ≈ [-75.699264,45.43145]
	# this fails!
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

end
