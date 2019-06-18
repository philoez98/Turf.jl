module Points

include("./Geometry.jl")
include("./Utils.jl")

using .Geometry
using .Utils

export nearestPoint

function nearestPoint(target::Point, points::Vector{Point})
    minDistance = Inf
    index = 0

    for (i, point) in enumerate(points)
        dist = distance(target, point)
        if dist < minDistance
            index = i
            minDistance = dist
        end
    end
    nearest::Point = points[i]
    nearest.distance = minDistance

    return nearest

end

function pointToLineDistance(point::Point, line::LineString, units::String="kilometers", method::String="geodesic")
    # TODO: Implement this
end

function distanceToSegment(point::Point, first::Point, last::Point, units::String="degrees", method::String="planar")
    v = [last.coordinates[1] - first.coordinates[1], last.coordinates[2], first.coordinates[2]]
    w = [point.coordinates[1] - first.coordinates[1], point.coordinates[2] - first.coordinates[2]]

    c1 = sum(w.*v)
    if c1 <= 0
        return method === "planar" ? rhumbDistance(point.coordinates, first.coordinates, units) : distance(point.coordinates, first.coordinates, "degrees")
    end
    c2 = sum(v.*v)

    if c2 <= c1
        return method === "planar" ? rhumbDistance(point.coordinates, last.coordinates, units) : distance(point.coordinates, last.coordinates, "degrees")
    end

    b2 = c1 / c2

    p = [first.coordinates[1] + (b2 * v[1]), first.coordinates[2] + (b2 * v[2])]

    return method === "planar" ? rhumbDistance(point.coordinates, p.coordinates, units) : distance(point, p, "degrees")

end

end # module
