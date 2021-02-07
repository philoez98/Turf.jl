@testset "booleans" begin
    line1 = LineString([[9.170356, 45.477985], [9.164434, 45.482551], [9.166644, 45.484003]])
    line2 = LineString([[9.169356, 45.477985], [9.163434, 45.482551], [9.165644, 45.484003]])

    line3 = LineString([[
						-111.544189453125,
						24.186847428521244
					],
					[
						-110.687255859375,
						24.966140159912975
					],
					[
						-110.4510498046875,
						24.467150664739002
					],
					[
						-109.9951171875,
						25.180087808990645
					]])

	line4 = LineString([[
						-111.4617919921875,
						24.05148034322011
					],
					[
						-110.8795166015625,
						24.681961205014595
					],
					[
						-110.841064453125,
						24.14174098050432
					],
					[
						-109.97863769531249,
						24.617057340809524
					]])

    @test parallel(line1, line2) == true
	@test parallel(line3, line4) == false

	poly1 = Polygon([[[0, 0], [1, 0], [1, 1], [0.5, 0.5], [0, 1], [0, 0]]])
	poly2 = Polygon([[[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]]])

	@test concave(poly1) == true
	@test concave(poly2) == false

	l1 = LineString([[0, 0], [1, 1], [1, 0], [0, 0]])
	l2 = LineString([[0, 0], [1, 0], [1, 1], [0, 0]])

	@test clockwise(l1) == true
	@test clockwise(l2) == false

	l3 = LineString([[0, 0], [3, 3], [4, 4]])
	p1 = Point([1,1])

	l4 = LineString([[0, 0], [3, 3]])
	p2 = Point([0, 0])

	p3 = Point([20, 20])
	l5 = LineString([[0, 0], [3, 3], [38.32, 5.96]])

	@test point_on_line(p2, l4, true) == false
	@test point_on_line(p3, l5, true) == false
	@test point_on_line(p1, l3) == true

	pt = Point([-77, 44])
	poly = Polygon([[ [-81, 41], [-81, 47], [-72, 47], [-72, 41], [-81, 41]]])

	@test point_in_polygon(pt, poly) == true

	poly3 = Polygon([[[1, 1], [1, 10], [10, 10], [10, 1], [1, 1]]])
	poly4 = Polygon([[[1, 1], [2, 2], [3, 2], [1, 1]]])
	line5 = LineString([[1, 1], [2, 3], [2, 3.5]])

	line6 = LineString([[1, 1], [1, 2], [1, 3], [1, 4]])
	poly5 = Polygon( [[[1, 1], [1, 20], [1, 3], [1, 4], [1, 1]]])
	line7 = LineString([[1, 2], [1, 3], [1, 3.5]])

	@test contains(poly3, poly4) == true
	@test contains(poly3, line5) == true
	@test contains(line6, Point([1, 2])) == true
	@test contains(poly3, poly5) == false
	@test contains(poly3 , line7) == false

	@test within(poly4, poly3) == true
	@test within(line5, poly3) == true
	@test within(poly5, poly3) == false
	@test within(Point([1, 2]), line6) == true
	@test within(line7, poly3) == false

	poly6 = Polygon([[[-11, -12], [-13, -12], [-13, -13], [-11, -13], [-11, -12]]])
	poly7 = Polygon([[[-1, 2], [3, 2], [3, 3], [-1, 3], [-1, 2]]])

	#@test disjoint(poly7, poly6) == true # <- fails
	@test disjoint(poly7, Point([1, 1])) == true
	@test disjoint(poly7, LineString([[0, 0], [12, 2], [12, 3], [12, 4]])) == true

	poly8 = Polygon([[[-1, 2], [-13, -12], [-13, -13], [-11, -13], [-1, 2]]])

	@test disjoint(poly8, poly7) == false

	line8 = LineString([[124.584961,-12.768946],[126.738281,-17.224758]])
	line9 = LineString([[123.354492,-15.961329],[127.22168,-14.008696]])

	@test line_intersects(line8, line9).coordinates ≈ [125.583754,-14.835723]

	line10 = LineString([
					[
						142.03125,
						-11.695273
					],
					[
						138.691406,
						-16.804541
					],
					[
						136.40625,
						-14.604847
					],
					[
						135.966797,
						-12.039321
					],
					[
						131.308594,
						-11.436955
					],
					[
						128.232422,
						-15.36895
					],
					[
						125.947266,
						-13.581921
					],
					[
						121.816406,
						-18.729502
					],
					[
						117.421875,
						-20.632784
					],
					[
						113.378906,
						-23.402765
					],
					[
						114.169922,
						-26.667096
					]])

	line11 = LineString([[
						117.861328,
						-15.029686
					],
					[
						122.124023,
						-24.886436
					],
					[
						132.583008,
						-22.309426
					],
					[
						132.890625,
						-7.754537]])

	points = line_intersects(line10, line11)
	@test points[1] ≈ [132.808697,-11.630938]
	@test points[2] ≈ [119.832884,-19.58857]

	@test crosses(LineString([[-2, 2], [4, 2]]), line6) == true
	@test crosses(LineString([[0.5, 2.5], [1, 1]]), poly7) == true
	@test crosses(MultiPoint([[1, 2], [12, 12]]), LineString([[1, 1], [1, 2], [1, 3], [1, 4]])) == true
	@test crosses(MultiPoint([[1, 0], [12, 12]]), LineString([[1, 1], [1, 2], [1, 3], [1, 4]])) == false
	@test crosses(LineString([[-2, 2], [-4, 2]]), poly7) == false


    pl1 = Polygon([[[0,0],[0,5],[5,5],[5,0],[0,0]]])
	pl2 = Polygon([[[1,1],[1,6],[6,6],[6,1],[1,1]]])

	@test overlap(pl1, pl2) == true
	@test_throws ErrorException overlap(pl1, Point([1, 1]))
	@test_throws ErrorException overlap(Point([1, 1]), pl2)

	pl3 = pl4 = Polygon([
					[
						[
							-53.57208251953125,
							28.287451910503744
						],
						[
							-53.33038330078125,
							28.29228897739706
						],
						[
							-53.34136962890625,
							28.430052892335723
						],
						[
							-53.57208251953125,
							28.287451910503744
						]
					]
				])
	@test overlap(pl3, pl4) == false

	mp1 = MultiPoint([
					[
						-36.05712890625,
						26.480407161007275
					],
					[
						-35.7220458984375,
						27.137368359795584
					],
					[
						-35.13427734375,
						26.83387451505858
					],
					[
						-35.4638671875,
						27.254629577800063
					],
					[
						-35.5462646484375,
						26.86328062676624
					],
					[
						-35.3924560546875,
						26.504988828743404
					]
				])
	mp2 = MultiPoint([
					[
						-35.4638671875,
						27.254629577800063
					],
					[
						-35.5462646484375,
						26.86328062676624
					],
					[
						-35.3924560546875,
						26.504988828743404
					],
					[
						-35.2001953125,
						26.12091815959972
					],
					[
						-34.9969482421875,
						26.455820238459893
					]
				])

	@test overlap(mp1, mp2) == true
	@test overlap(mp1, mp2) == overlap(mp2, mp1)
end
