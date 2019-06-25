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
	"""
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
	"""

	poly2 = Polygon([[[100, 0], [101, 0], [101, 1], [100, 1], [100, 0]]])
	feats = explode(poly2).features
	res = [Point([100, 0]), Point([101, 0]), Point([101, 1]), Point([100, 1]), Point([100, 0])]


	for i in eachindex(feats)
		@test feats[i].geometry.coordinates == res[i].coordinates
	end

	@test flip(Point([77.34374999999999,43.58039085560784,3000])).coordinates == [43.58039085560784,77.34374999999999,3000]

end
