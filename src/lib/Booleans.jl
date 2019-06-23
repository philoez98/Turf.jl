"""Takes a ring and return true or false whether or not the ring is clockwise or counter-clockwise."""
function clockwise(line::Union{LineString, Vector{Position}})::Bool

    let ring end

    if geotype(line) === :LineString
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
function concave(poly::Polygon)
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


function equal(geo1::T, geo2::T) where {T <: AbstractGeometry}
    geotype(geo1) !== geotype(geo2) && return false

    geotype(geo1) === :Point && return comparePoints(geo1.coordinates, geo2.coordinates)
    geotype(geo1) === :LineString && return compareLines(geo1.coordinates, geo2.coordinates)
end

function comparePoints(p1::Position, p2::Position)
    length(p1) !== length(p2) && return false

    for i in eachindex(p1)
        round(p1[i]; digits=10) !== round(p2[i]; digits=10) && return false
    end

    return true
end

function compareLines(p1::Vector{Position}, p2::Vector{Position})
    # TODO: complete this
    length(p1[1]) !== length(p2[1]) && return false
end

"""Return `true` if each segment of `line1` is parallel to the correspondent segment of `line2`"""
function parallel(line1::LineString, line2::LineString)::Bool
    seg1 = lineSegment(line1)
    seg2 = lineSegment(line2)

    for i in eachindex(seg1)
        coors2 = nothing
        coors1 = seg1[i].coordinates

        try
            coors2 = seg2[i].coordinates
        catch e
            isa(e, BoundsError) && break

        end
        isParallel(coors1, coors2) == false && return false
    end

    return true
end

"""Compare slopes"""
@inline function isParallel(p1::Vector{Position}, p2::Vector{Position})
    slope1 = bearingToAzimuth(rhumbBearing(p1[1], p2[1]))
    slope2 = bearingToAzimuth(rhumbBearing(p1[2], p2[2]))

    return slope1 === slope2
end
