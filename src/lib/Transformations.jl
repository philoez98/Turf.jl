const origin_options = ["sw", "se", "nw", "ne", "center", "centroid"]

"""
Rotates any geojson Feature or Geometry of a specified angle, around its `centroid` or a given `pivot` point;
all rotations follow the right-hand rule.
"""
function transformRotate(; geojson::T, angle::Real, pivot::Point=nothing, mutate::Bool=false) where {T <: AbstractGeometry}
    if angle === 0
        return geojson
    end

    geo = []

    if mutate
        geo = deepcopy(geojson)
    else
        geo = geojson
    end

    type = geotype(geo)

    pivot === nothing && (pivot = centroid(geo))
    coords = geo.coordinates

    if type === :Point

        initAngle = rhumbBearing(pivot.coordinates, coords)
        finalAngle = initAngle + angle
        dist = rhumbDistance(pivot.coordinates, coords)
        newCoords = rhumbDestination(pivot.coordinates, dist, finalAngle).coordinates
        coords[1] = newCoords[1]
        coords[2] = newCoords[2]

    elseif type === :Polygon || type === :MultiLineString

        for i in eachindex(coords[1])
            initAngle = rhumbBearing(pivot.coordinates, coords[1][i])
            finalAngle = initAngle + angle
            dist = rhumbDistance(pivot.coordinates, coords[1][i])
            newCoords = rhumbDestination(pivot.coordinates, dist, finalAngle).coordinates
            coords[1][i][1] = newCoords[1]
            coords[1][i][2] = newCoords[2]
        end
    elseif type === :LineString

        for i in eachindex(coords)
            initAngle = rhumbBearing(pivot.coordinates, coords[i])
            finalAngle = initAngle + angle
            dist = rhumbDistance(pivot.coordinates, coords[i])
            newCoords = rhumbDestination(pivot.coordinates, dist, finalAngle).coordinates
            coords[i][1] = newCoords[1]
            coords[i][2] = newCoords[2]
        end
    end

    return geo
end



"""
Scale a GeoJson from a given point by a factor of scaling (ex: factor=2 would make the GeoJson 200% larger).
If a FeatureCollection is provided, the origin point will be calculated based on each individual Feature.
"""
function transformScale(geojson::T, factor::Float64, origin::String="centroid") where {T <: AbstractFeatureCollection}

    for (i, feat) in enumerate(geojson.features)
        geojson.features[i] = scale(feat, factor, origin)
    end

    return geojson
end


"""
Scale a Feature.
"""
function scale(feature::Feature, factor::Real, origin::String="centroid")

    if !(origin in origin_options)
        throw(error("'$(origin)' is not a valid option. Allowed values are: $(origin_options)"))
    end

    (geotype(feature.geometry) === :Point || factor === 1) && return feature

    type = geotype(feature.geometry)
    coords = feature.geometry.coordinates
    center = getOrigin(feature, origin)

    if type === :LineString

        for i in eachindex(coords)
            start = rhumbDistance(center.coordinates, coords[i])
            bearing = rhumbBearing(center.coordinates, coords[i])
            distance = start * factor

            newCoords = rhumbDestination(center.coordinates, distance, bearing).coordinates

            coords[i][1] = newCoords[1]
            coords[i][2] = newCoords[2]

            if length(coords) === 3
                coords[i][3] *= factor
            end
        end
    elseif type === :Polygon || type === :MultiLineString

        for i in eachindex(coords[1])
            start = rhumbDistance(center.coordinates, coords[1][i])
            bearing = rhumbBearing(center.coordinates, coords[1][i])
            distance = start * factor

            newCoords = rhumbDestination(center.coordinates, distance, bearing).coordinates

            coords[1][i][1] = newCoords[1]
            coords[1][i][2] = newCoords[2]

            if length(coords) === 3
                coords[1][i][3] *= factor
            end
        end
    end

    return feature

end

function getOrigin(geojson::Feature, origin::String)
    if origin == nothing || origin === ""
        origin = "centroid"
    end

    sw = ["sw", "southwest", "westsouth", "bottomleft"]
    se = ["se", "southeast", "eastsouth", "bottomright"]
    nw = ["nw", "northwest", "westnorth", "topleft"]
    ne = ["ne", "northeast", "eastnorth", "topright"]

    box = bbox(geojson)
    west = box[1]
    south = box[2]
    east = box[3]
    north = box[4]

    if origin in sw
        return Point([west, south])
    elseif origin in se
        return Point([east, south])
    elseif origin in nw
        return Point([west, north])
    elseif origin in ne
        return Point([east, north])
    elseif origin === "center"
        return center(geojson.geometry)
    elseif origin === "centroid"
        return centroid(geojson.geometry)
    else
        throw(error("Invalid origin."))
    end
end

"""
Takes a Geometry or a FeatureCollection and returns all positions as Points.
"""
function explode(geojson::T, pointsOnly::Bool=false)::FeatureCollection where {T <: Union{AbstractFeatureCollection, AbstractGeometry}}
    points::Vector{Point} = []

    if geotype(geojson) === :FeatureCollection
        for feat in geojson.features
            geom = feat.geometry

            if geotype(geom) === :Point
                push!(points, geom)

            elseif geotype(geom) === :Polygon || geotype(geom) === :MultiLineString

                for i in eachindex(geom.coordinates[1])
                    push!(points, Point(geom.coordinates[1][i]))
                end
            elseif geotype(geom) === :LineString

                for i in eachindex(geom.coordinates)
                    push!(points, Point(geom.coordinates[i]))
                end

            end
        end
    else
        if geotype(geojson) === :Point
            push!(points, geojson)

        elseif geotype(geojson) === :Polygon || geotype(geojson) === :MultiLineString

            for i in eachindex(geojson.coordinates[1])
                push!(points, Point(geojson.coordinates[1][i]))
            end
        elseif geotype(geojson) === :LineString

            for i in eachindex(geojson.coordinates)
                push!(points, Point(geojson.coordinates[i]))
            end

        end

    end

    pointsOnly && return points
    return FeatureCollection([Feature(x) for x in points])
end

"""Take input Features and Geometries and flips all of their coordinates from `[x, y]` to `[y, x]`."""
function flip(geojson::T, mutate::Bool=false) where {T <: Union{AbstractFeature, AbstractGeometry}}
    type = geotype(geojson)
    coords = []
    geom = nothing

    !mutate && (geojson = deepcopy(geojson))

    if type === :Feature
        geom = geojson.geom
        coords = geojson.geometry.coordinates
    else
        geom = geojson
        coords = geojson.coordinates
    end

    if geotype(geom) === :Point
        x = coords[1]
        y = coords[2]
        coords[1] = y
        coords[2] = x
    elseif geotype(geom) === :Polygon || geotype(geom) === :MultiLineString
        for i in eachindex(coords)
            x = coords[i][1]
            y = coords[i][2]
            coords[i][1] = y
            coords[i][2] = x
        end
    elseif geotype(geom) === :LineString
        for i in eachindex(coords[1])
            x = coords[1][i][1]
            y = coords[1][i][2]
            coords[1][i][1] = y
            coords[1][i][2] = x
        end
    end

    return geojson

end

"""
Cohen-Sutherland line clippign algorithm, adapted to efficiently to
handle polylines rather than just segments
"""
function lineclip(points::Vector{P}, bbox::Vector{T}) where {T <: Real, P <: AbstractPoint}
    len = length(points)
    codeA = bitcode(points[1], bbox)
    parts = []
    results = []
    let a, b, codeB, lastCode end

    for i in 2:len
        a = points[i - 1]
        b = points[i]

        codeB = lastCode = bitcode(b, bbox)

        while true
            if !(codeA | codeB)
                push!(parts, a)

                if codeB !== lastCode
                    push!(parts, b)

                    if i < len - 1
                        push!(results, parts)
                        parts = []
                    end
                elseif i === len - 1
                    push!(parts, b)
                end

                break
            elseif (codeA & codeB)
                break

            elseif codeA
                a = intersection(a, b, codeB, bbox)
                codeA = bitcode(b, bbox)

            else
                b = intersection(a, b, codeB, bbox)
                codeB = bitcode(b, bbox)
            end
        end
        codeA = lastCode
    end
    length(part) > 0 && push!(results, parts)

    return results
end

"""Sutherland-Hodgeman polygon clipping algorithm."""
function polygonclip(points::Vector{P}, bbox::Vector{T}) where {P <: AbstractPoint, T <: Real}
    let result, edge, prev, pInside, p, inside end
    ranges = [1, 2, 4, 8]

    for edge in ranges
        result = []

        prev = points[length(points) - 1]
        pInside = !(bitcode(prev, bbox) & edge)

        for i in 1:length(points)
            p = points[i]
            inside = !(bitcode(p, bbox) & edge)

            inside !== pInside && push!(result, intersection(prev, p, edge, bbox))

            inside && push!(result, p)

            prev = p
            pInside = inside
        end

        points = result
        length(points) <= 0 && break
    end

    return result
end


"""Intersect a segment against one of the 4 lines that make up the bbox"""
function intersection(a::Point, b::Point, edge, bbox::Vector{T}) where {T <: Real}
    return edge & 8 ? [a[1] + (b[1] - a[1]) * (bbox[4] - a[1]) / (b[1] - a[1]), bbox[4]] : # top
           edge & 4 ? [a[1] + (b[1] - a[1]) * (bbox[2] - a[2]) / (b[2] - a[2]), bbox[2]] : # bottom
           edge & 2 ? [bbox[3], a[2] + (b[2] - a[2]) * (bbox[3] - a[1]) / (b[1] - a[1])] : # right
           edge & 1 ? [bbox[1], a[2] + (b[2] - a[2]) * (bbox[1] - a[1]) / (b[1] - a[1])] : # left
           nothing
end


@inline function bitcode(p::Point, bbox::Vector{T}) where {T <: Real}
    code = 0

    p[1] < bbox[1] && (code |= 1)
    p[1] > bbox[3] && (code |= 2)

    p[2] < bbox[2] && (code |= 4)
    p[2] > bbox[4] && (code |= 8)

    return code
end

""" Finds the tangents of a Polygon from a Point."""
function polygonTangents(pt::Point, poly::Polygon)
    ptCoords = pt.coordinates
    polyCoords = poly.coordinates

    box = bbox(poly)
    nearest = nothing
    nearestIndex = 1

    if ptCoords[1] > box[1] && ptCoords[1] < box[3] && ptCoords[2] > box[2] && ptCoords[2] < box[4]
        nearest = nearestPoint(pt, explode(poly, true))
        nearestIndex = findfirst(x -> x == nearest, explode(poly, true))
    end

    rtan = polyCoords[1][nearestIndex]
    ltan = polyCoords[1][1]

    if nearest != nothing
        nearest.coordinates[2] < ptCoords[2] && (ltan = polyCoords[1][nearestIndex])
    end

    enext = nothing
    eprev = isleft(polyCoords[1][1], polyCoords[1][length(polyCoords[1]) - 1], ptCoords)
    out = process(polyCoords[1], ptCoords, eprev, enext, rtan, ltan)

    rtan = out[1]
    ltan = out[2]

    return FeatureCollection([Feature(Point(rtan)), Feature(Point(ltan))])
end

function isleft(point1, point2, point3)
    return (point2[1] - point1[1]) * (point3[2] - point1[2]) - (point3[1] - point1[1]) * (point2[2] - point1[2])
end

function isabove(point1, point2, point3)
    return isleft(point1, point2, point3) > 0
end

function isbelow(point1, point2, point3)
    return isleft(point1, point2, point3) < 0
end

function process(polyCoords, ptCoords, eprev, enext, rtan, ltan)
    for i in 1:length(polyCoords) - 1
        curr = polyCoords[i]
        next = polyCoords[i + 1]

        (i === length(polyCoords) - 1) && (next = polyCoords[1])

        enext = isleft(curr, next, ptCoords)
        if eprev <= 0 && enext > 0
            !(isbelow(ptCoords, curr, rtan)) && (rtan = curr)
        elseif eprev > 0 && enext <= 0
            !(isabove(ptCoords, curr, ltan)) && (ltan = curr)
        end
        eprev = enext
    end
    return [rtan, ltan]
end
