@inline function createSegments(coords::Vector{Position})
    segments = []
    if length(coords) === 2
        push!(segments, Linestring([coords[1], coords[2]]))
        return segments
    end

    for i in 1:length(coords)-1
        line = LineString([coords[i], coords[i + 1]])
        push!(segments, line)
    end
    return segments
end

function linesegmentFeature(geojson::Union{LineString, Polygon}, result::Vector{LineString})
    coords = geojson.coordinates

    segments = createSegments(coords)

    for s in segments
        push!(result, s)
    end

    return result
end

"""
    linesegment(geojson::LineString)

Create a 2-vertex LineString segments from a LineString.

# Examples
```jldoctest
julia> line = LineString([[0, 1], [2, 3], [3, 3]])
LineString(Array{Float64,1}[[0.0, 1.0], [2.0, 3.0], [3.0, 3.0]])

julia> linesegment(line)
2-element Array{LineString,1}:
 LineString(Array{Float64,1}[[0.0, 1.0], [2.0, 3.0]])
 LineString(Array{Float64,1}[[2.0, 3.0], [3.0, 3.0]])
```
"""
function linesegment(geojson::LineString)
    result::Vector{LineString} = []
    linesegmentFeature(geojson, result)

    return result
end

function to360(a::Real)
    b = a % 360
    b < 0 && (b += 360)

    return b
end

"""
    linearc(center::Point, radius::Real, bearing1::Real, bearing2::Real, steps::Real=64., units::String="kilometers")

Create a circular arc, of a circle of the given radius and center point, between bearing1 and bearing2;
0 bearing is North of center point, positive clockwise.

# Examples
```jldoctest
julia> point = Point([35, 45])
Point([35.0, 45.0])

julia> linearc(point, 5, 0, 10, 5)
LineString(Array{Float64,1}[[35.0, 45.045], [35.0111, 45.0443]])

julia> linearc(point, 5, 0, 270, 5)
LineString(Array{Float64,1}[[35.0, 45.045], [35.0605, 45.0139], [35.0374, 44.9636], [34.9626, 44.9636], [34.9364, 45.0]])
```
"""
function linearc(center::Point, radius::Real, bearing1::Real, bearing2::Real, steps::Real=64., units::String="kilometers")
    angle1 = to360(bearing1)
    angle2 = to360(bearing2)

    angle1 === angle2 && return LineString(circle(center=center, radius=radius, steps=steps, units=units).coordinates[1])

    startdeg = angle1
    enddeg = (angle1 < angle2) ?  angle2 :  angle2 + 360
    α = startdeg

    coords = []
    i = 0

    while α < enddeg
        push!(coords, destination(center.coordinates, radius, α).coordinates)
        i += 1
        α = startdeg + i * 360 / steps
    end

    α > enddeg && push!(coords, destination(center.coordinates, radius, enddeg).coordinates)

    return LineString(coords)
end
