"""Takes a ring and return true or false whether or not the ring is clockwise or counter-clockwise."""
function clockwise(line::Union{LineString, Vector{Position}})::Bool

    let ring end
    # geotype fails with Vector{Vector{...}}
    if geotype(line) === :LineString
        ring = line.coordinates
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
    coords = poly.coordinates

    length(coords[1]) <= 4 && return false

    sign = false
    n = length(coords[1]) - 1

    for i in 1:n
        j = ((i + 1) % n) === 0 ? 1 : (i + 1) % n
        m = ((i + 2) % n) === 0 ? 1 : (i + 2) % n

        dx1 = coords[1][m][1] - coords[1][j][1]
        dy1 = coords[1][m][2] - coords[1][j][2]
        dx2 = coords[1][i][1] - coords[1][j][1]
        dy2 = coords[1][i][2] - coords[1][j][2]

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

"""
    pointOnLine(point::Point, line::LineString, ignoreEndVertices::Bool=false)::Bool

Returns true if a point is on a line. Accepts a optional parameter to ignore the
start and end vertices of the linestring.
"""
function pointOnLine(point::Point, line::LineString, ignoreEndVertices::Bool=false)::Bool
    pCoords = point.coordinates
    lCoords = line.coordinates

    ignore = "none"
    for i in 1:length(lCoords) - 1
        if ignoreEndVertices == true
            i === 1 && (ignore = "start")
            i === length(lCoords) - 2 && (ignore = "end")
            (i === 1 && i + 1 === length(lCoords) - 1) && (ignore = "both")
        end
        isPointOnSegment(lCoords[i], lCoords[i + 1], pCoords, ignore) && return true
    end
    return false
end

@inline function isPointOnSegment(start::Position, stop::Position, coords::Position, excludeBoundary::String="none")::Bool
    x, y = coords
    x1, y1 = start
    x2, y2 = stop

    dxc = x - x1
    dyc = y - y1
    dx1 = x2 - x1
    dy1 = y2 - y1

    cross = dxc * dy1 - dyc * dx1
    cross != 0 && return false

    if excludeBoundary === "none"
        if abs(dx1) >= abs(dy1)
            return dx1 > 0 ? x1 <= x && x <= x2 : x2 <= x && x <= x1
        end
        return dy1 > 0 ? y1 <= y && y <= y2 : y2 <= y && y <= y1
    elseif excludeBoundary === "start"
        if abs(dx1) >= abs(dy1)
             return dx1 > 0 ? x1 < x && x <= x2 : x2 <= x && x < x1
        end
        return dy1 > 0 ? y1 < y && y <= y2 : y2 <= y && y < y1
    elseif excludeBoundary === "end"
        if abs(dx1) >= abs(dy1)
            return dx1 > 0 ? x1 <= x && x < x2 : x2 < x && x <= x1
        end
        return dy1 > 0 ? y1 <= y && y < y2 : y2 < y && y <= y1
    elseif excludeBoundary === "both"
        if abs(dx1) >= abs(dy1)
            return dx1 > 0 ? x1 < x && x < x2 : x2 < x && x < x1
        end
        return dy1 > 0 ? y1 < y && y < y2 : y2 < y && y < y1
    end
    return false
end


"""
Takes a Point and a Polygon and determines if the point
resides inside the polygon. The polygon can be convex or concave. The function accounts for holes.
"""
function pointInPolygon(point::Point, polygon::Union{Polygon, MultiPolygon}, ignoreBoundary::Bool=false)

    pt = point.coordinates
    poly = polygon.coordinates

    inBBox(pt, bbox(polygon)) == false && return false

    geotype(polygon) === :Polygon && (poly = [poly])

    inside = false
    for i in eachindex(poly)
        if inRing(pt, poly[i][1], ignoreBoundary)
            inHole = false
            k = 1

            while k < length(poly[i]) && !inHole
                inRing(pt, poly[i][k], !ignoreBoundary) == true && (inHole = true)
                k += 1
            end

            !inHole && (inside = true)
        end
    end
    return inside
end

function inRing(pt::Position, ring::Vector{Position}, ignoreBoundary::Bool=false)
    inside = false

    (ring[1][1] == ring[length(ring) - 1][1] && ring[1][2] == ring[length(ring) - 1][1]) && (ring = ring[1, length(ring) - 1])

    for i in 1:length(ring) - 1
        j = i + 1

        xi = ring[i][1]
        yi = ring[i][2]
        xj = ring[j][1]
        yj = ring[j][2]

        onBoundary = (pt[2] * (xi - xj) + yi * (xj - pt[1]) + yj * (pt[1] - xi) == 0) &&
            ((xi - pt[1]) * (xj - pt[1]) <= 0) && ((yi - pt[2]) * (yj - pt[2]) <= 0)

        onBoundary && return !ignoreBoundary

        intersect = ((yi > pt[2]) !== (yj > pt[2])) && (pt[1] < (xj - xi) * (pt[2] - yi) / (yj - yi) + xi)

        intersect && (inside =  !inside)
    end

    return inside
end

function inBBox(pt::Position, bbox::Vector{Float64})
    return bbox[1] <= pt[1] &&  bbox[2] <= pt[2] &&
        bbox[3] >= pt[1] && bbox[4] >= pt[2]
end

"""
Return True if the second geometry is completely contained by the first geometry.
The interiors of both geometries must intersect and, the interior and boundary of the secondary (geometry b)
must not intersect the exterior of the primary (geometry a).
`contains` returns the exact opposite result of `within`.
"""
function contains(ft1::AbstractGeometry, ft2::AbstractGeometry)::Bool
    type1 = geotype(ft1)
    type2 = geotype(ft2)

    coords1 = ft1.coordinates
    coords2 = ft2.coordinates

    if type1 === :Point
        if type2 === :Point
            return coords1[1] == coords2[1] && coords1[2] == coords2[2]
        end
        throw(error("$(type2) is not a supported type."))
    elseif type1 === :LineString
        if type2 === :Point
            return pointOnLine(ft2, ft1, true)
        elseif type2 === :LineString
            return lineOnLine(ft1, ft2)
        else
            throw(error("$(type2) is not a supported type."))
        end

    elseif type1 === :Polygon
        if type2 === :Point
            return pointInPolygon(ft2, ft1, true)
        elseif type2 === :LineString
            return lineInPolygon(ft1, ft2)
        elseif type2 === :Polygon
            return polygonInPolygon(ft1, ft2)
        else
            throw(error("$(type2) is not a supported type."))
        end
    else
        throw(error("$(type1) is not a supported type."))
    end
end

contains(ft1::Feature, ft2::Feature) = contains(ft1.geometry, ft2.geometry)

function lineInPolygon(poly::Polygon, line::LineString)
    out = false

    polybox = bbox(poly)
    linebox = bbox(line)

    !(bboxOverlap(polybox, linebox)) && return false

    coords = line.coordinates

    for i in 1:length(coords) - 1
        mid = [(coords[i][1] + coords[i + 1][1]) / 2, (coords[i][2] + coords[i + 1][2]) / 2]
        if pointInPolygon(Point(mid), poly, true)
            out = true
            break
        end
    end
    return out
end

function lineOnLine(line1::LineString, line2::LineString)
    for i in eachindex(line1.coordinates)
        !(pointOnLine(line1.coordinates[i], line2)) && return false
    end
    return true
end

function polygonInPolygon(ft1::Polygon, ft2::Polygon, reverse::Bool=false)
    polybox1 = bbox(ft1)
    polybox2 = bbox(ft2)
    coords = []

    if reverse
        coords = ft1.coordinates
        !(bboxOverlap(polybox2, polybox1)) && return false

        for ring in coords
            for coord in ring
                !(pointInPolygon(Point(coord), ft2)) && return false
            end
        end
    else
        coords = ft2.coordinates
        !(bboxOverlap(polybox1, polybox2)) && return false

        for ring in coords
            for coord in ring
                !(pointInPolygon(Point(coord), ft1)) && return false
            end
        end
    end

    return true
end

function bboxOverlap(box1::Vector{T}, box2::Vector{T}) where {T <: Real}
    box1[1] > box2[1] && return false
    box1[3] < box2[3] && return false
    box1[2] > box2[2] && return false
    box1[4] < box2[4] && return false
    return true
end


function within(ft1::AbstractGeometry, ft2::AbstractGeometry)::Bool
    type1 = geotype(ft1)
    type2 = geotype(ft2)

    coords1 = ft1.coordinates
    coords2 = ft2.coordinates

    if type1 === :Point
        if type2 === :LineString
            return pointOnLine(ft1, ft2, true)
        elseif type2 === :Polygon
            return pointInPolygon(ft1, ft2, true)
        end
        throw(error("$(type2) is not a supported type."))
    elseif type1 === :LineString
        if type2 === :Polygon
            return lineInPolygon(ft1, ft2)
        elseif type2 === :LineString
            return lineOnLine(ft1, ft2)
        else
            throw(error("$(type2) is not a supported type."))
        end

    elseif type1 === :Polygon
        if type2 === :Polygon
            return polygonInPolygon(ft1, ft2, true)
        else
            throw(error("$(type2) is not a supported type."))
        end
    else
        throw(error("$(type1) is not a supported type."))
    end
end

function lineInPolygon(line::LineString, poly::Polygon)
    polybox = bbox(poly)
    linebox = bbox(line)

    !(bboxOverlap(polybox, linebox)) && return false

    coords = line.coordinates
    inside = false

    for i in 1:length(coords) - 1
        !(pointInPolygon(Point(coords[i]), poly)) && return false
        !inside && (inside = pointInPolygon(Point(coords[i]), poly, true))
        if !inside
            mid = [(coords[i][1] + coords[i + 1][1]) / 2, (coords[i][2] + coords[i + 1][2]) / 2]
            inside = pointInPolygon(Point(mid), poly, true)
        end
    end
    return inside
end

"""Return `true` if the intersection of the two geometries is an empty set."""
function disjoint(geom1::AbstractGeometry, geom2::AbstractGeometry)
    type1 = geotype(geom1)
    type2 = geotype(geom2)

    coords1 = geom1.coordinates
    coords2 = geom2.coordinates

    if type1 === :Point
        if type2 === :Point
            return !(coords1[1] == coords2[1] && coords1[2] == coords2[2])
        elseif type2 === :LineString
            return !pointOnLine(geom1, geom2)
        elseif type2 === :Polygon
            !pointInPolygon(geom1, geom2)
        end

    elseif type1 === :LineString
        if type2 === :Point
            return !pointOnLine(geom2, geom1)
        elseif type2 === :LineString
            return !lineOnLine(geom1, geom2)
        elseif type2 === :Polygon
            return !lineInPolygon(geom2, geom1)
        end
    elseif type1 === :Polygon
        if type2 === :Point
            return !pointInPolygon(geom2, geom1)
        elseif type2 === :LineString
            return !lineInPolygon(geom1, geom2)
        elseif type2 === :Polygon
            return !polyInPoly(geom2, geom1)
        end
    end
    return false
end

disjoint(geom1::AbstractFeature, geom2::AbstractFeature) = disjoint(geom1.geometry, geom2.geometry)

function polyInPoly(poly1::Polygon, poly2::Polygon)
    coords1 = poly1.coordinates
    coords2 = poly2.coordinates

    for ring in coords1
        for coord in ring
            (pointInPolygon(Point(coord), poly2)) && return true
        end
    end

    for ring in coords2
        for coord in ring
            (pointInPolygon(Point(coord), poly1)) && return true
        end
    end

    inter = lineIntersects(polygonToLine(poly1), polygonToLine(poly2))
    inter != nothing && return true

    return false

end

"""Find a point that intersects LineStrings with two coordinates each."""
function lineIntersects(line1::AbstractLineString, line2::AbstractLineString)
    coords1 = line1.coordinates
    coords2 = line2.coordinates

    (length(coords1) == 2 && length(coords2) == 2) && return intersects(line1, line2)

    result = []

    for i in 1:length(coords1) - 1
        for j in 1:length(coords2) - 1
            inter = intersects(LineString([coords1[i], coords1[i + 1]]), LineString([coords2[j], coords2[j + 1]]))
            inter != nothing && push!(result, inter.coordinates)
        end
    end
    return unique!(result)
end




function intersects(line1::AbstractLineString, line2::AbstractLineString)
    coords1 = line1.coordinates
    coords2 = line2.coordinates

    (length(coords1) != 2 || length(coords2) != 2) && throw(error("Lines must contain only 2 coordinates."))

    x1, y1 = coords1[1]
    x2, y2 = coords1[2]
    x3, y3 = coords2[1]
    x4, y4 = coords2[2]

    d = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1))
    a = ((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3))
    b = ((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3))

    if d == 0
        if a == 0 && b == 0
            return nothing
        end
        return nothing
    end

    ã = a / d
    b̃ = b / d

    if ã >= 0 && ã <= 1 && b̃ >= 0 && b̃ <= 1
        x = x1 + (ã * (x2 - x1))
        y = y1 + (ã * (y2 - y1))
        return Point([x, y])
    end

    return nothing
end

"""
Return `true` if the intersection results in a geometry whose dimension is one less than
the maximum dimension of the two source geometries and the intersection set is interior to
both source geometries.
"""
function crosses(ft1::AbstractGeometry, ft2::AbstractGeometry)::Bool
    type1 = geotype(ft1)
    type2 = geotype(ft2)

    if type1 === :MultiPoint
        if type2 === :LineString
            return MpCrossLs(ft1, ft2)
        elseif type2 === :Polygon
            return MpCrossPoly(ft1, ft2)
        else
            throw(error("Geometry $(type2) is not supported."))
        end
    elseif type1 === :LineString
        if type2 === :MultiPoint
            return MpCrossLs(ft2, ft1)
        elseif type2 === :Polygon
            return LsCrossPoly(ft1, ft2)
        elseif type2 === :LineString
            return LsCrossLs(ft1, ft2)
        else
            throw(error("Geometry $(type2) is not supported."))
        end
    elseif type1 === :Polygon
        if type2 === :MultiPoint
            return MpCrossPoly(ft2, ft1)
        elseif type2 === :LineString
            return LsCrossPoly(ft2, ft1)
        else
            throw(error("Geometry $(type2) is not supported."))
        end
    end

    throw(error("Geometry $(type1) is not supported."))
end

crosses(ft1::Feature, ft2::Feature) = crosses(ft1.geometry, ft2.geometry)


function MpCrossLs(geom1::MultiPoint, geom2::LineString)
    intPoint = false
    extPoint = false

    pLength = length(geom1.coordinates)
    i = 1

    while i < pLength && !intPoint && !extPoint

        for j in 1:length(geom2.coordinates) - 1
            incVertices = "both"

            (j === 1 || j === length(geom2.coordinates) - 2) && (incVertices = "none")

            if isPointOnSegment(geom2.coordinates[j], geom2.coordinates[j + 1], geom1.coordinates[i], incVertices)
                intPoint = true
            else
                extPoint = true

            end

        end
        i += 1
    end

    return intPoint && extPoint
end

function LsCrossLs(line1::LineString, line2::LineString)
    inter = lineIntersects(line1, line2)

    if length(inter) > 0
        for i in 1:length(line1.coordinates) - 1
            for j in 1:length(line2.coordinates) - 1
                incVertices = "both"

                (j === 1 || j === length(line2.coordinates) - 2) && (incVertices = "none")

                isPointOnSegment(line1.coordinates[i], line1.coordinates[i + 1], line2.coordinates[j], incVertices) && return true
            end
        end
    end
    return false
end

function LsCrossPoly(line::LineString, poly::Polygon)
    line2 = polygonToLine(poly)
    inter = lineIntersects(line, line2)

    (length(inter) > 0) && return true

    return false
end

function MpCrossPoly(mp::MultiPoint, poly::Polygon)
    intPoint = false
    extPoint = false
    pLength = length(mp.coordinates[1])
    i = 1

    while i < pLength && intPoint && extPoint
        if pointInPolygon(Point(mp.coordinates[1][i]), poly)
            intPoint = true
        else
            extPoint = true
        end
        i += 1
    end
    return intPoint && extPoint
end
