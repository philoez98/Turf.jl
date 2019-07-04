"""
Create a 2-vertex LineString segments from a LineString.
"""
function linesegment(geojson::LineString)
    result::Vector{LineString} = []
    linesegmentFeature(geojson, result)

    return result
end

function linesegmentFeature(geojson::Union{LineString, Polygon}, result::Vector{LineString})
    coords = geojson.coordinates

    segments = createSegments(coords)

    for s in segments
        push!(result, s)
    end

    return result
end

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

"""
    linearc(center::Point, radius::Real, bearing1::Real, bearing2::Real, steps::Real=64., units::String="kilometers")
    
Creates a circular arc, of a circle of the given radius and center point, between bearing1 and bearing2;
0 bearing is North of center point, positive clockwise.
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

function to360(a::Real)
    b = a % 360
    b < 0 && (b += 360)

    return b
end
