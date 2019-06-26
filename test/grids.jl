@testset "grids" begin
    bbox1 = [-96.6357421875,31.12819929911196,-84.9462890625,40.58058466412764]
    bbox2 = [-220.78125,-80.64703474739618,-29.53125,78.34941069014629]
    bbox3 = [-5.1,-5.1,5.1,5.1]

    feat2 = pointGrid(bbox2, 500., nothing, "miles")
    @test feat2.features[1].geometry.coordinates ≈ [-199.504452, -77.132893]

    feat3 = rectangleGrid(bbox3, 1., 1., nothing, "degrees")
    @test feat3.features[1].geometry.coordinates ≈ [
					[
						[
							-5.025791,
							-5.005842
						],
						[
							-5.025791,
							-4.004674
						],
						[
							-4.020633,
							-4.004674
						],
						[
							-4.020633,
							-5.005842
						],
						[
							-5.025791,
							-5.005842
						]
					]] atol=0.0001

    feat1 = hexGrid(bbox1, 25., nothing, false, "miles")
"""
	@test feat1.features[1].geometry.coordinates ≈
    [
                [
                    [
                        -95.704593,
                        31.467449
                    ],
                    [
                        -95.927937,
                        31.780802
                    ],
                    [
                        -96.374626,
                        31.780802
                    ],
                    [
                        -96.59797,
                        31.467449
                    ],
                    [
                        -96.374626,
                        31.154096
                    ],
                    [
                        -95.927937,
                        31.154096
                    ],
                    [
                        -95.704593,
                        31.467449
                    ]
                ]]
				"""
end
