using Random: shuffle, RandomDevice

"""
    circle(; center::Union{Point, Position}, radius::Real=5., steps::Integer=64, units::String="kilometers")

Take a Point or a Position and calculate the circle polygon given a radius in degrees, radians, miles, or kilometers; and steps for precision.

# Examples
```jldoctest
julia> point = Point([35, 45])
Point([35.0, 45.0])

julia> circle(center=point, steps=5)
Polygon(Array{Array{Float64,1},1}[[[34.9395, 45.0139], [34.9626, 44.9636], [35.0374, 44.9636], [35.0605, 45.0139], [35.0, 45.045], [34.9395, 45.0139]]])

julia> circle(center=point, radius=2.5, steps=5, units="degrees")
Polygon(Array{Array{Float64,1},1}[[[31.5893, 45.7231], [32.99, 42.9571], [37.01, 42.9571], [38.4107, 45.7231], [35.0, 47.5029], [31.5893, 45.7231]]])
```
"""
function circle(; center::Union{Point, Position}, radius::Real=5., steps::Integer=64, units::String="kilometers")
    coords = []
    position::Position = []

    if geotype(center) === :Point
        position = Position(center.coordinates)
    else
        position = center
    end

    for i in 1:steps
        push!(coords, destination(position, Float64(radius), i * (-360 / steps), units).coordinates)
    end
    push!(coords, coords[1])

    return Polygon([coords])
end

"""
    sector(center::Point, radius::Real, bearing1::Real, bearing2::Real, steps::Integer=64, units::String="kilometers")

Creates a circular sector of a circle of given radius and center Point,
between (clockwise) bearing1 and bearing2; 0 bearing is North of center point, positive clockwise.

# Examples
```jldoctest
julia> point = Point([35, 45])
Point([35.0, 45.0])

julia> sector(point, 5, 0, 0.5, 5)
Polygon(Array{Array{Float64,1},1}[[[35.0, 45.0], [35.0, 45.045], [35.0006, 45.045], [35.0, 45.0]]])

julia> sector(point, 5, 0, 270, 5)
Polygon(Array{Array{Float64,1},1}[[[35.0, 45.0], [35.0, 45.045], [35.0605, 45.0139], [35.0374, 44.9636], [34.9626, 44.9636], [34.9364, 45.0], [35.0, 45.0]]])
```
"""
function sector(center::Point, radius::Real, bearing1::Real, bearing2::Real, steps::Integer=64, units::String="kilometers")
    complem(bearing1) === complem(bearing2) && return circle(center=center, radius=radius, steps=steps, units=units)

    coords = center.coordinates
    arc = linearc(center, radius, bearing1, bearing2, steps, units)
    sliceCoords = [[coords]]

    for i in eachindex(arc.coordinates)
        push!(sliceCoords[1], arc.coordinates[i])
    end
    push!(sliceCoords[1], coords)

    return Polygon(sliceCoords)
end

function complem(a::Real)
    b = a % 360
    b < 0 && (b += 360)

    return b
end


"""
    smallest_circle(points::Vector{Point})::Polygon

Create the minimum bounding circle for the given sets of points recursively using [Welzl's algorithm](https://en.wikipedia.org/wiki/Smallest-circle_problem)
"""
function smallest_circle(points::Vector{Point})::Polygon
    random = shuffle(RandomDevice(), points)
    r = -1
    c = circle(center=Point([0, 0]), radius=r)

    for i in 1:length(random)-1
        p = random[i]
        println(p)
        (r < 0 || !point_in_polygon(p, c)) && (c = one_point_circle(points, i+1, p))
    end
    return c
end



function circumcircle(a::Position, b::Position, c::Position)
    ox = min(min(a[1], b[1]), c[1]) + max(min(a[1], b[1]), c[1]) / 2
    oy = min(min(a[2], b[2]), c[2]) + max(min(a[2], b[2]), c[2]) / 2

    ax = a[1] - ox
    ay = a[2] - oy
    bx = b[1] - ox
    by = b[2] - oy
    cx = c[1] - ox
    cy = c[2] - oy

    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    d == 0 && return circle(center=Point([0, 0]), radius=-1), -1, Point([0, 0])

    x = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d
    y = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d

    p = Point([x + ox, y + oy])
    da = hypot(p.coordinates[1] - a[1], p.coordinates[2] - a[2])
    db = hypot(p.coordinates[1] - b[1], p.coordinates[2] - b[2])
    dc = hypot(p.coordinates[1] - c[1], p.coordinates[2] - c[2])
    r = max(max(da, db), dc)

    return circle(center=p, radius=r, steps=100), r, p
end

function one_point_circle(points::Vector{Point}, size::Integer, center::Point)
    r = 1
    c = circle(center=center, radius=r)

    for i in 1:size
        p = points[i]
        if point_in_polygon(p, c)
            if r == 0
                x = (center.coordinates[1] + p.coordinates[1]) / 2
                y = (center.coordinates[2] + p.coordinates[2]) / 2

                q = Point([x, y])
                println("p: $(p), q: $(q), cen: $(center)")
                r = diameter(p, center, q)
                c = circle(center=q, radius=r, steps=100)
            else
                c = two_point_circle(points, i, center, p)
            end
        end
    end

    return c
end

function two_point_circle(points::Vector{Point}, size::Integer, center::Point, p::Point)
    x = (center.coordinates[1] + p.coordinates[1]) / 2
    y = (center.coordinates[2] + p.coordinates[2]) / 2
    println("cen : $center, p: $p")
    q = Point([x, y])
    r = diameter(p, q, center)
    circ = circle(center=q, radius=r, steps=100)

    rl = rr = -1
    left_c = right_c = Point([0, 0])

    left = circle(center=left_c, radius=rl)
    right = circle(center=right_c, radius=rr)

    d_pq = subtract(q, p)

    for i in 1:size
        t = points[i]

        point_in_polygon(t, circ) && continue

        d = cross(d_pq, subtract(t, center))
        c, r, new_center = circumcircle(center.coordinates, p.coordinates, t.coordinates)

        if r < 0
            continue
        elseif d > 0 && ( rl < 0 || cross(d_pq, subtract(new_center, center)) > cross(d_pq, subtract(left_c, center)))
            left = c
            rl = r
        elseif d < 0 || ( rr < 0 || cross(d_pq, subtract(new_center, center)) > cross(d_pq, subtract(right_c, center)))
            right = c
            rr = r
        end
    end

    (rl < 0 && rr < 0) && return circ
    rl < 0 && return right
    rr < 0 && return left

    return rl <= rr ? left : right
end



function diameter(a::Point, b::Point, c::Point)
    da = hypot(c.coordinates[1] - a.coordinates[1], c.coordinates[2] - a.coordinates[2])
    db = hypot(c.coordinates[1] - b.coordinates[1], c.coordinates[2] - b.coordinates[2])
    println("MAX: $(max(da, db))")
    return max(da, db)
end

function subtract(p::Point, q::Point)
    return Point([p.coordinates[1] - q.coordinates[1],
        p.coordinates[2] - q.coordinates[2]])
end

function cross(p::Point, q::Point)
    return p.coordinates[1] * q.coordinates[2] -
        p.coordinates[2] * q.coordinates[1]
end
