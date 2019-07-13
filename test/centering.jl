@testset "centering" begin

    p1 = Point([4.831961989402771, 45.75764678012361])
    l1 = LineString([[
        4.86020565032959,
        45.76884015325622
      ],
      [
        4.85994815826416,
        45.749558161214516
      ]])

    poly = Polygon([[[
          4.8250579833984375,
          45.79398056386735
        ],
        [
          4.882392883300781,
          45.79254427435898
        ],
        [
          4.910373687744141,
          45.76081677972451
        ],
        [
          4.894924163818359,
          45.7271539426975
        ],
        [
          4.824199676513671,
          45.71337148333104
        ],
        [
          4.773387908935547,
          45.74021417890731
        ],
        [
          4.778022766113281,
          45.778418789239055
        ],
        [
          4.8250579833984375,
          45.79398056386735
        ]]])

    @test centroid(p1).coordinates == [4.831961989402771, 45.75764678012361]
    @test centroid(l1).coordinates == [4.860076904296875,45.75919915723537]
    @test centroid(poly).coordinates == [4.839177131652832,45.76256007199914]

    poly2 = Polygon([[[0, 0], [0, 0], [0, 0], [0, 0]]])

    @test meancenter(poly).coordinates ≈ [4.839177,45.76256]
    @test meancenter(poly2).coordinates == [0, 0]
    @test meancenter(l1).coordinates ≈ [4.860077,45.759199]

    @test masscenter(poly2).coordinates == [0, 0]
    @test masscenter(poly).coordinates == [4.840728965137111,45.75581209996416]

    coll = FeatureCollection([Feature(Point([0, 0])),
      Feature(Point([1, 0])), Feature(Point([0, 1])),
      Feature(Point([5, 8]))])

    coll2 = FeatureCollection([Feature(Point([0, 0])),Feature(Point([9, 9])), Feature(Point([9.25, 9.25])),
      Feature(Point([9.5, 9.5])), Feature(Point([9.75, 9.75])), Feature(Point([10, 10]))])

    @test mediancenter(coll).coordinates ≈ [0.383876,0.616989] atol=0.0001
    @test mediancenter(coll2).coordinates ≈ [9.254246,9.254246] atol=0.0001

    fc = FeatureCollection([Feature(Point([4.833351373672485,45.760809294695534])),
      Feature(Point([4.8331475257873535,45.760296567821456])), Feature(Point([4.833984374999999,45.76073818687033])),
      Feature(Point([4.834005832672119,45.76022171678877]))])

    @test centroid(fc).coordinates ≈ [4.8336222767829895,45.76051644154402] atol=0.001
end
