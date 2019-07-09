"""
    bezier(line::LineString; resolution::Real=10000, sharpness::Real=0.85)

Take a LineString and returns a curved version
by applying a [Bezier spline](http://en.wikipedia.org/wiki/B%C3%A9zier_spline)
algorithm.
"""
function bezier(line::LineString; resolution::Real=10000, sharpness::Real=0.85)

    points = []

    for p in line.coordinates
        push!(points, p)
    end

    spline = Splines.Spline(points=points, duration=resolution, sharpness=sharpness)

    coords = []
    for i in 1:10:spline.duration
        pos = Splines.pos(i)

        floor(i / 100) % 2 == 0 && push!(coords, [pos.coordinates[1], pos.coordinates[2]])

    end

    return LineString(coords)

end
