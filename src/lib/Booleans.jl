include("../geojson/Geometries.jl")


"""Takes a ring and return true or false whether or not the ring is clockwise or counter-clockwise."""
function clockwise(line::Union{Geometries.LineString, Vector{Geometries.Position}})::Bool

    let ring end

    if typeof(line) <: Geometries.AbstractLineString
        ring = line.geometry.coordinates
    else
        ring = line[1]
    end

    sum = 0
    i = 1

    prev = []
    cur = []
    while i < length(ring)
        prev = isempty(cur) ? ring[1] : cur
        cur = ring[i]
        sum += (cur[1] - prev[1]) * (cur[2] + prev[2])
        i += 1
    end

    return sum > 0
end

"""Takes a polygon and return true or false as to whether it is concave or not."""
function concave(ploy::Geometries.Polygon)
    coords = poly.geometry.coordinates

    length(coords[1]) <= 4 && return false

    sign = false
    n = length(coords[1])

    for i in 1:n
        dx1 = coords[1][(i + 2) % n][1] - coords[1][(i + 1) % n][1]
        dy1 = coords[1][(i + 2) % n][2] - coords[1][(i + 1) % n][2]
        dx2 = coords[1][i][1] - coords[1][(i + 1) % n][1]
        dy2 = coords[1][i][2] - coords[1][(i + 1) % n][2]

        cross = (dx1 * dy2) - (dy1 * dx2)

        if i === 0
            sign = cross > 0
        elseif sign !== (cross > 0)
            return true
        end
    end

    return false
end


function equal(geo1::Geometries.Geometry, geo2::Geometries.Geometry)
    geo1.type !== geo2.type && return false

    geo1.type === "Point" && return comparePoints(geo1.coordinates, geo2.coordinates)
    geo1.type === "LineString" && return compareLines(geo1.coordinates, geo2.coordinates)
end

function comparePoints(p1::Geometries.Position, p2::Geometries.Position)
    length(p1) !== length(p2) && return false

    for i in eachindex(p1)
        round(p1[i]; digits=10) !== round(p2[i]; digits=10) && return false
    end

    return true
end

function compareLines(p1::Vector{Geometries.Position}, p2::Vector{Geometries.Position})
    # TODO: complete this
    length(p1[1]) !== length(p2[1]) && return false
end
