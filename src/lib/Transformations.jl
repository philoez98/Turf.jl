const origin_options = ["sw", "se", "nw", "ne", "center", "centroid"]

"""
    transform_rotate([geojson::T[, angle::Real], pivot::Point=nothing, mutate::Bool=false) where {T <: AbstractGeometry}

Rotates any geojson Feature or Geometry of a specified angle, around its `centroid` or a given `pivot` point;
all rotations follow the right-hand rule.

# Examples
```jldoctest
julia> point = Point([-75.69926351308823,45.43145021122502])
Point([-75.6993, 45.4315])

julia> transform_rotate(geojson=point, angle=80., pivot=Point([-75.6, 45.3]))
Point([-75.433, 45.3915])
```
"""
function transform_rotate(; geojson::T, angle::Real, pivot::Point=nothing, mutate::Bool=false) where {T <: AbstractGeometry}
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

        initAngle = rhumb_bearing(pivot.coordinates, coords)
        finalAngle = initAngle + angle
        dist = rhumb_distance(pivot.coordinates, coords)
        newCoords = rhumb_destination(pivot.coordinates, dist, finalAngle).coordinates
        coords[1] = newCoords[1]
        coords[2] = newCoords[2]

    elseif type === :Polygon || type === :MultiLineString

        for i in eachindex(coords)
            for j in eachindex(coords[i])
                initAngle = rhumb_bearing(pivot.coordinates, coords[i][j])
                finalAngle = initAngle + angle
                dist = rhumb_distance(pivot.coordinates, coords[i][j])
                newCoords = rhumb_destination(pivot.coordinates, dist, finalAngle).coordinates
                coords[i][j][1] = newCoords[1]
                coords[i][j][2] = newCoords[2]
            end
        end
    elseif type === :LineString

        for i in eachindex(coords)
            initAngle = rhumb_bearing(pivot.coordinates, coords[i])
            finalAngle = initAngle + angle
            dist = rhumb_distance(pivot.coordinates, coords[i])
            newCoords = rhumb_destination(pivot.coordinates, dist, finalAngle).coordinates
            coords[i][1] = newCoords[1]
            coords[i][2] = newCoords[2]
        end
    end

    return geo
end

transform_rotate!(; geojson::T, angle::Real, pivot::Point=nothing) where {T <: AbstractGeometry} = transform_rotate(geojson=geojson, angle=angle, pivot=pivot, mutate=true)

"""
    transform_translate([geojson::T[, distance::R[, direction::R]]], vertical::R=0, mutate::Bool=false, units::String="kilometers") where {T <: Union{AbstractFeature, AbstractGeometry}, R <: Real}

Moves any geojson Feature or Geometry of a specified distance along a Rhumb Line
on the provided direction angle.

# Examples
```jldoctest
julia> poly = Polygon([[[0, 29], [3.5, 29], [2.5, 32], [0, 29]]])
Polygon(Array{Array{Float64,1},1}[[[0.0, 29.0], [3.5, 29.0], [2.5, 32.0], [0.0, 29.0]]])

julia> transform_translate(poly, 300, 70)
Polygon(Array{Array{Float64,1},1}[[[2.91184, 29.9228], [6.41184, 29.9228], [5.50479, 32.9228], [2.91184, 29.9228]]])
```
"""
function transform_translate(geojson::T, distance::R, direction::R, vertical::R=0, mutate::Bool=false, units::String="kilometers") where {T <: Union{AbstractFeature, AbstractGeometry}, R <: Real}
    (distance == 0 && vertical == 0) && return geojson

    if distance < 0
        distance  = -distance
        direction = -direction
    end

    !mutate && (geojson = deepcopy(geojson))

    coords = []
    type = nothing

    if typeof(geojson) <: AbstractGeometry
        coords = geojson.coordinates
        type = geotype(geojson)
    else
        coords = geojson.geometry.coordinates
        type = geotype(geojson.geometry)
    end

    if type === :Point
        newCoords = rhumb_destination(coords, distance, direction, units).coordinates
        coords[1] = newCoords[1]
        coords[2] = newCoords[2]
        (vertical != 0 && length(coords) === 3) && (coords[3] += vertical)

    elseif type === :Polygon || type === :MultiLineString
        for i in eachindex(coords)
            for j in eachindex(coords[i])
                newCoords = rhumb_destination(coords[i][j], distance, direction, units).coordinates
                coords[i][j][1] = newCoords[1]
                coords[i][j][2] = newCoords[2]
                (vertical != 0 && length(coords[i][j]) === 3) && (coords[i][j][3] += vertical)
            end
        end
    elseif type === :LineString
        for i in eachindex(coords)
            newCoords = rhumb_destination(coords[i], distance, direction, units).coordinates
            coords[i][1] = newCoords[1]
            coords[i][2] = newCoords[2]
            (vertical != 0 && length(coords[i]) === 3) && (coords[i][3] += vertical)
        end
    end

    return geojson
end

transform_translate!(geojson::T, distance::R, direction::R, vertical::R=0, units::String="kilometers") where {T <: Union{AbstractFeature, AbstractGeometry}, R <: Real} = transform_translate(geojson, distance, direction, vertical, true, units)


"""
    transform_scale([geojson::T[, factor::Float64]], origin::String="centroid") where {T <: AbstractFeatureCollection}

Scale a GeoJson from a given point by a factor of scaling (ex: factor=2 would make the GeoJson 200% larger).
If a FeatureCollection is provided, the origin point will be calculated based on each individual Feature.

# Examples
```jldoctest
julia> coll = FeatureCollection([Feature(Point([-75.69926351308823,45.43145021122502])), Feature(Polygon([[[0, 29], [3.5, 29], [2.5, 32], [0, 29]]]))])
FeatureCollection{Feature}(Feature[Feature(Point([-75.6993, 45.4315]), Dict{String,Any}()), Feature(Polygon(Array{Array{Float64,1},1}[[[0.0, 29.0], [3.5, 29.0], [2.5, 32.0], [0.0, 29.0]]]), Dict{String,Any}())], nothing, nothing)

julia> transform_scale(coll, 0.1)
FeatureCollection{Feature}(Feature[Feature(Point([-75.6993, 45.4315]), Dict{String,Any}()), Feature(Polygon(Array{Array{Float64,1},1}[[[1.3495, 29.675], [1.70067, 29.675], [1.59896, 29.975], [1.3495, 29.675]]]), Dict{String,Any}())], nothing, nothing)
```
"""
function transform_scale(geojson::T, factor::Real, origin::String="centroid", mutate::Bool=false) where {T <: AbstractFeatureCollection}

    geo = []
    if !mutate
        geo = deepcopy(geojson)
    else
        geo = geojson
    end

    for (i, feat) in enumerate(geo.features)
        geo.features[i] = scale(feat, factor, origin)
    end

    return geo
end

transform_scale!(geojson::T, factor::Real, origin::String="centroid") where {T <: AbstractFeatureCollection} = transform_scale(geojson, factor, origin, true)

"""
    scale([feature::Feature[, factor::Real]], origin::String="centroid")

Scale a Feature.

# Examples
```jldoctest
julia> feature = Feature(Polygon([[[0, 29], [3.5, 29], [2.5, 32], [0, 29]]]))
Feature(Polygon(Array{Array{Float64,1},1}[[[0.0, 29.0], [3.5, 29.0], [2.5, 32.0], [0.0, 29.0]]]), Dict{String,Any}())

julia> scale(feature, 0.1)
Feature(Polygon(Array{Array{Float64,1},1}[[[1.3495, 29.675], [1.70067, 29.675], [1.59896, 29.975], [1.3495, 29.675]]]), Dict{String,Any}())
```
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
            start = rhumb_distance(center.coordinates, coords[i])
            bearing = rhumb_bearing(center.coordinates, coords[i])
            distance = start * factor

            newCoords = rhumb_destination(center.coordinates, distance, bearing).coordinates

            coords[i][1] = newCoords[1]
            coords[i][2] = newCoords[2]

            if length(coords) === 3
                coords[i][3] *= factor
            end
        end
    elseif type === :Polygon || type === :MultiLineString

        for i in eachindex(coords)
            for j in eachindex(coords[i])
                start = rhumb_distance(center.coordinates, coords[i][j])
                bearing = rhumb_bearing(center.coordinates, coords[i][j])
                distance = start * factor

                newCoords = rhumb_destination(center.coordinates, distance, bearing).coordinates

                coords[i][j][1] = newCoords[1]
                coords[i][j][2] = newCoords[2]

                if length(coords) === 3
                    coords[i][j][3] *= factor
                end
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
    explode([geojson::T], pointsOnly::Bool=false)::FeatureCollection where {T <: Union{AbstractFeatureCollection, AbstractGeometry}}

Takes a Geometry or a FeatureCollection and returns all positions as Points.

# Examples
```jldoctest
julia> poly = Polygon([[[100, 0], [101, 0], [101, 1], [100, 1], [100, 0]]])
Polygon(Array{Array{Float64,1},1}[[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]]])

julia> explode(poly, true)
5-element Array{Point,1}:
 Point([100.0, 0.0])
 Point([101.0, 0.0])
 Point([101.0, 1.0])
 Point([100.0, 1.0])
 Point([100.0, 0.0])

```
"""
function explode(geojson::T, pointsOnly::Bool=false) where {T <: Union{AbstractFeatureCollection, AbstractGeometry}}
    points::Vector{Point} = []

    if geotype(geojson) === :FeatureCollection
        for feat in geojson.features
            geom = feat.geometry

            if geotype(geom) === :Point
                push!(points, geom)

            elseif geotype(geom) === :Polygon || geotype(geom) === :MultiLineString

                for i in eachindex(geom.coordinates)
                    for j in eachindex(geom.coordinates[i])
                        push!(points, Point(geom.coordinates[i][j]))
                    end
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

            for i in eachindex(geojson.coordinates)
                for j in eachindex(geojson.coordinates[i])
                    push!(points, Point(geojson.coordinates[i][j]))
                end
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

"""
    flip([geojson::T], mutate::Bool=false) where {T <: Union{AbstractFeature, AbstractGeometry}}

Take input Features and Geometries and flips all of their coordinates from `[x, y]` to `[y, x]`.

# Examples
```jldoctest
julia> point = Point([77.34374999999999,43.58039085560784,3000])
Point([77.3437, 43.5804, 3000.0])

julia> flip(point)
Point([43.5804, 77.3437, 3000.0])
```
"""
function flip(geojson::T, mutate::Bool=false) where {T <: Union{AbstractFeature, AbstractGeometry}}
    type = geotype(geojson)
    coords = []
    geom = nothing

    geo = []
    if !mutate
        geo = deepcopy(geojson)
    else
        geo = geojson
    end

    if type === :Feature
        geom = geo.geom
        coords = geo.geometry.coordinates
    else
        geom = geo
        coords = geo.coordinates
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

    return geo

end

flip!(geojson::T) where {T <: Union{AbstractFeature, AbstractGeometry}} = flip(geojson, true)

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

                if codeB != lastCode
                    push!(parts, b)

                    if i < len
                        push!(results, parts)
                        parts = []
                    end
                elseif i === len
                    push!(parts, b)
                end
                break
            elseif (codeA & codeB)
                break

            elseif codeA
                a = intersection(a, b, codeB, bbox)
                codeA = bitcode(a, bbox)

            else
                b = intersection(a, b, codeB, bbox)
                codeB = bitcode(b, bbox)
            end
        end
        codeA = lastCode
    end
    length(parts) > 0 && push!(results, parts)

    return results
end

"""
Sutherland-Hodgeman polygon clipping algorithm.
"""
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


"""
Intersect a segment against one of the 4 lines that make up the bbox
"""
function intersection(a::Point, b::Point, edge, bbox::Vector{T}) where {T <: Real}
    ac = a.coordinates
    bc = b.coordinates
    return edge & 8 ? [ac[1] + (bc[1] - ac[1]) * (bbox[4] - ac[2]) / (bc[2] - ac[2]), bbox[4]] : # top
           edge & 4 ? [ac[1] + (bc[1] - ac[1]) * (bbox[2] - ac[2]) / (bc[2] - ac[2]), bbox[2]] : # bottom
           edge & 2 ? [bbox[3], ac[2] + (bc[2] - ac[2]) * (bbox[3] - ac[1]) / (bc[1] - ac[1])] : # right
           edge & 1 ? [bbox[1], ac[2] + (bc[2] - ac[2]) * (bbox[1] - ac[1]) / (bc[1] - ac[1])] : # left
           nothing
end


@inline function bitcode(p::Point, bbox::Vector{T}) where {T <: Real}
    code = 0
    c = p.coordinates

    c[1] < bbox[1] && (code |= 1)
    c[1] > bbox[3] && (code |= 2)

    c[2] < bbox[2] && (code |= 4)
    c[2] > bbox[4] && (code |= 8)

    return code
end

"""
    polygon_tangents(pt::Point, poly::Polygon)

Finds the tangents of a Polygon from a Point.

# Examples
```jldoctest
julia> point = Point([92.46093749999999,54.67383096593114])
Point([92.4609, 54.6738])

julia> poly = Polygon([[[100, 0], [101, 0], [101, 1], [100, 1], [100, 0]]])
Polygon(Array{Array{Float64,1},1}[[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]]])

julia> poly = Polygon([[[48.1641, 20.6328], [76.6406, 20.6328], [76.6406, 38.8226], [48.1641, 38.8226], [48.1641, 20.6328]]])
Polygon(Array{Array{Float64,1},1}[[[48.1641, 20.6328], [76.6406, 20.6328], [76.6406, 38.8226], [48.1641, 38.8226], [48.1641, 20.6328]]])

julia> polygon_tangents(point, poly)
FeatureCollection{Feature}(Feature[Feature(Point([48.1641, 38.8226]), Dict{String,Any}()), Feature(Point([76.6406, 20.6328]), Dict{String,Any}())], nothing, nothing)
 ```
"""
function polygon_tangents(pt::Point, poly::Polygon)
    ptCoords = pt.coordinates
    polyCoords = poly.coordinates

    box = bbox(poly)
    nearest = nothing
    nearestIndex = 1

    if ptCoords[1] > box[1] && ptCoords[1] < box[3] && ptCoords[2] > box[2] && ptCoords[2] < box[4]
        nearest = nearestpoint(pt, explode(poly, true))
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


"""
    polygon_to_line([poly::Polygon])

Converts a Polygon to LineString or MultiLineString

# Examples
```jldoctest
julia> poly = Polygon([[[-2.275543, 53.464547],[-2.275543, 53.489271],[-2.215118, 53.489271],[-2.215118, 53.464547],[-2.275543, 53.464547]]])
Polygon(Array{Array{Float64,1},1}[[[-2.27554, 53.4645], [-2.27554, 53.4893], [-2.21512, 53.4893], [-2.21512, 53.4645], [-2.27554, 53.4645]]])

julia> polygon_to_line(poly)
LineString(Array{Float64,1}[[-2.27554, 53.4645], [-2.27554, 53.4893], [-2.21512, 53.4893], [-2.21512, 53.4645], [-2.27554, 53.4645]])
```
"""
function polygon_to_line(poly::Polygon)
    return coordinatesToLine(poly.coordinates)
end

function coordinatesToLine(coords::Vector{Vector{Position}})
    length(coords) > 1 && return MultiLineString(coords)
    return LineString(coords[1])
end

"""
    convert_to([geojson::AbstractGeometry[, projection::String]], mutate::Bool=false)

Convert a GeoJSON geojsonect to the defined `projection`.
"""
function convert_to(geojson::AbstractGeometry, projection::String, mutate::Bool=false)
    allowed_proj = ["mercator", "wgs84"]

    !(projection in allowed_proj) && throw(error("$(projection) is not a valid option."))

    type = geotype(geojson)

    type === :Point && return Point((projection === "mercator") ? to_mercator(geojson) : to_WGS84(geojson))

    geo = []
    if !mutate
        geo = deepcopy(geojson)
    else
        geo = geojson
    end

    coords = geo.coordinates
    if type === :LineString
        for i in eachindex(coords)
            newCoords = (projection === "mercator") ? to_mercator(Point(coords[i])) : to_WGS84(Point(coords[i]))
            coords[i][1] = newCoords[1]
            coords[i][2] = newCoords[2]
        end
    elseif type === :Polygon
        for i in eachindex(coords[1])
            newCoords = (projection === "mercator") ? to_mercator(Point(coords[1][i])) : to_WGS84(Point(coords[1][i]))
            coords[1][i][1] = newCoords[1]
            coords[1][i][2] = newCoords[2]
        end
    else
        throw(error("$(type) is not supported."))
    end

    return geo
end

convert_to(geojson::Feature, projection::String, mutate::Bool=false) = convert_to(geojson.geometry, projection, mutate)
convert_to!(geojson::Feature, projection::String) = convert_to(geojson.geometry, projection, true)
convert_to!(geojson::AbstractGeometry, projection::String) = convert_to(geojson, projection, true)

"""
    combine([ft::FeatureCollection])

Combine a FeatureCollection of Point, LineString, or Polygon features
into MultiPoint, MultiLineString, or MultiPolygon features.

# Examples
```jldoctest
julia> l1 = LineString([[102.0,-10.0],[130.0,4.0]])
LineString(Array{Float64,1}[[102.0, -10.0], [130.0, 4.0]])

julia> l2 = LineString([[40.0,-20.0],[150.0,18.0]])
LineString(Array{Float64,1}[[40.0, -20.0], [150.0, 18.0]])

julia> combine(FeatureCollection([Feature(l1), Feature(l2)]))
FeatureCollection{Feature}(Feature[Feature(MultiLineString(Array{Array{Float64,1},1}[[[102.0, -10.0], [130.0, 4.0], [40.0, -20.0], [150.0, 18.0]]]), Dict{String,Any}())], nothing, nothing)
```
"""
function combine(ft::FeatureCollection)

    groups = Dict{String, Any}([
        "MultiPoint" => [[], Dict{String, Any}()],
        "MultiLineString" => [[], Dict{String, Any}()],
        "MultiPolygon" => [[], Dict{String, Any}()]]
    )

    function add_to_group(ft::Feature, key::String, multi::Bool)
        if multi
            if geotype(ft.geometry) === :MultiPoint
                for i in eachindex(ft.geometry.coordinates)
                    push!(groups[key][1], ft.geometry.coordinates[i])
                end
            elseif geotype(ft.geometry) === :MultiLineString
                for i in eachindex(ft.geometry.coordinates[1])
                    push!(groups[key][1], ft.geometry.coordinates[1][i])
                end
            else
                for i in eachindex(ft.geometry.coordinates[1][1])
                    push!(groups[key][1], ft.geometry.coordinates[1][1][i])
                end
            end
        else
            if geotype(ft.geometry) === :Point
                groups[key][1] = vcat(groups[key][1], [ft.geometry.coordinates])
            elseif geotype(ft.geometry) === :LineString
                groups[key][1] = vcat(groups[key][1], ft.geometry.coordinates)
            elseif geotype(ft.geometry) === :Polygon
                groups[key][1] = vcat(groups[key][1], ft.geometry.coordinates...)
            end
        end # function

        merge!(groups[key][2], ft.properties)
    end

    for feat in ft.features
        type = geotype(feat.geometry)
        try
            add_to_group(feat, string(type), true)
        catch e
            if isa(e, KeyError)
                newType = "Multi" * string(type)
                add_to_group(feat, newType, false)
            else
                throw(e)
            end
        end
    end

    features::Vector{Feature} = []

    groups = filter( x -> length(groups[x.first][1]) > 0, groups)

    for g in keys(groups)
        if g === "MultiPoint"
            point = MultiPoint(values(groups[g][1]))
            push!(features, Feature(point, groups[g][2]))
        elseif g === "MultiLineString"
            point = MultiLineString([values(groups[g][1])])
            push!(features, Feature(point, groups[g][2]))
        elseif g === "MultiPolygon"
            point = MultiPolygon([[values(groups[g][1])]])
            push!(features, Feature(point, groups[g][2]))
        end
    end

    return FeatureCollection(features)
end

"""
    tag(fc1::FeatureCollection, fc2::FeatureCollection, in_field::String, out_field::String)::FeatureCollection

Take a set of Points and a set of Polygons and performs a spatial join.
"""
function tag(fc1::FeatureCollection, fc2::FeatureCollection, in_field::String, out_field::String)::FeatureCollection
    points = deepcopy(fc1)
    polys = deepcopy(fc2)

    for feat in points.features
        geotype(feat.geometry) !== :Point && throw(error("$(feat.geometry) is not a supported geometry."))
        for feat2 in polys.features
            geotype(feat2.geometry) !== :Polygon && throw(error("$(feat2.geometry) is not a supported geometry."))

            if point_in_polygon(feat.geometry, feat2.geometry)
                feat.properties[out_field] = feat2.properties[in_field]
            end
        end
    end
    return points
end


"""
    simplify(geojson::FeatureCollection, tolerance::Real=1., hq::Bool=false, mutate::Bool=false)

Take a FeatureCollection and return a simplified version.
Internally uses an adaptation of [simplify-js](http://mourner.github.io/simplify-js/) to perform simplification using the Ramer-Douglas-Peucker algorithm.

# Examples
```jldoctest
julia> poly = Polygon([[
    [-70.603637, -33.399918],
    [-70.614624, -33.395332],
    [-70.639343, -33.392466],
    [-70.659942, -33.394759],
    [-70.683975, -33.404504],
    [-70.697021, -33.419406],
    [-70.701141, -33.434306],
    [-70.700454, -33.446339],
    [-70.694274, -33.458369],
    [-70.682601, -33.465816],
    [-70.668869, -33.472117],
    [-70.646209, -33.473835],
    [-70.624923, -33.472117],
    [-70.609817, -33.468107],
    [-70.595397, -33.458369],
    [-70.587158, -33.442901],
    [-70.587158, -33.426283],
    [-70.590591, -33.414248],
    [-70.594711, -33.406224],
    [-70.603637, -33.399918]]])

julia> simplify(poly, 0.01)
Polygon(Array{Array{Float64,1},1}[[[-70.6036, -33.3999], [-70.684, -33.4045], [-70.7011, -33.4343], [-70.6943, -33.4584], [-70.6689, -33.4721], [-70.6098, -33.4681], [-70.5872, -33.4429], [-70.6036, -33.3999]]])
```
"""
function simplify(geojson::FeatureCollection, tolerance::Real=1., hq::Bool=false, mutate::Bool=false)
    tolerance < 0 && throw(error("Invalid tolerance. It must be a positive number!"))

    !mutate && (geojson = deepcopy(geojson))

    for feat in geojson.features
        geometry = feat.geometry

        simplify_geometry(geometry, tolerance, hq)
    end

    return geojson
end


"""
    simplify(geojson::FeatureCollection, tolerance::Real=1., hq::Bool=false)

Take a FeatureCollection and return a simplified version, modifying the original FeatureCollection.
Internally uses an adaptation of [simplify-js](http://mourner.github.io/simplify-js/) to perform simplification using the Ramer-Douglas-Peucker algorithm.
"""
simplify!(geojson::FeatureCollection, tolerance::Real=1., hq::Bool=false) = simplify(geojson, tolerance, hq, true)


"""
    simplify(geojson::AbstractGeometry, tolerance::Real=1., hq::Bool=false, mutate::Bool=false)

Take a GeoJSON Geometry and return a simplified version.
Internally uses an adaptation of [simplify-js](http://mourner.github.io/simplify-js/) to perform simplification using the Ramer-Douglas-Peucker algorithm.
"""
function simplify(geojson::AbstractGeometry, tolerance::Real=1., hq::Bool=false, mutate::Bool=false)
    !mutate && (geojson = deepcopy(geojson))
    return simplify_geometry(geojson, tolerance, hq)
end


"""
    simplify(geojson::AbstractGeometry, tolerance::Real=1., hq::Bool=false)

Take a GeoJSON Geometry and return a simplified version, modifying the original Geometry.
Internally uses an adaptation of [simplify-js](http://mourner.github.io/simplify-js/) to perform simplification using the Ramer-Douglas-Peucker algorithm.
"""
simplify!(geojson::AbstractGeometry, tolerance::Real=1., hq::Bool=false) = simplify(geojson, tolerance, hq, true)


"""
    simplify(geojson::AbstractFeature, tolerance::Real=1., hq::Bool=false, mutate::Bool=false)

Take a GeoJSON Feature and return a simplified version.
Internally uses an adaptation of [simplify-js](http://mourner.github.io/simplify-js/) to perform simplification using the Ramer-Douglas-Peucker algorithm.
"""
simplify(geojson::AbstractFeature, tolerance::Real=1., hq::Bool=false, mutate::Bool=false) = simplify(geojson.geometry, tolerance, hq, mutate)


"""
    simplify(geojson::AbstractFeature, tolerance::Real=1., hq::Bool=false)

Take a GeoJSON Feature and return a simplified version, modifying the original Feature.
Internally uses an adaptation of [simplify-js](http://mourner.github.io/simplify-js/) to perform simplification using the Ramer-Douglas-Peucker algorithm.
"""
simplify!(geojson::AbstractFeature, tolerance::Real=1., hq::Bool=false) = simplify(geojson, tolerance, hq, true)


function simplify_geometry(geom::AbstractGeometry, tol::Real, hq::Bool)
    type = geotype(geom)

    isequal(type, :Point) || isequal(type, :MultiPoint) && return geom

    clean!(geom)

    coords = geom.coordinates

    if isequal(type, :LineString)
        geom.coordinates = simplify_line(coords, tol, hq)
    elseif isequal(type, :MultiLineString)
        geom.coordinates = map(lines -> simplify_line(lines, tol, hq), coords[1])

    elseif isequal(type, :Polygon)
        geom.coordinates = simplify_polygon(coords, tol, hq)
    else
        geom.coordinates = map(rings -> simplify_polygon(rings, tol, hq), coords[1])
    end

    return geom
end


function simplify_line(coords, tol::Real, hq::Bool)
    return map(coords -> (coords[3]) ? [coords[1], coords[2], coords[3]] : [coords[1], coords[2]],
        DP_simplify(map(coord -> [coord[1], coord[2], coord[3]], coords), tol, hq))
end

function simplify_polygon(coords, tol::Real, hq::Bool)
    return map(ring -> begin
        pts = map(coord -> [coord[1], coord[2]], ring)

        length(pts) < 4 && throw(error("Invalid Polygon."))

        simple = map(p -> [p.coordinates[1], p.coordinates[2]], DP_simplify([Point(p) for p in pts], tol, hq))

        while !valid(simple)
            tol -= tol * 0.01
            simple = map(p -> [p.coordinates[1], p.coordinates[2]], DP_simplify([Point(p) for p in pts], tol, hq))
        end

        ((simple[length(simple)][1] != simple[1][1]) ||
            (simple[length(simple)][2] != simple[1][2])) && push!(simple, simple[1])

        return simple

        end, coords)
end

function valid(ring)
    length(ring) < 3 && return false

    return !(isequal(length(ring), 3) && (isequal(ring[3][1], ring[1][1]) && isequal(ring[3][2], ring[1][2])))
end

function dist²(p1::Point, p2::Point)
    dx = p1.coordinates[1] - p2.coordinates[1]
    dy = p1.coordinates[2] - p2.coordinates[2]

    return dx * dx + dy * dy
end

function segdist²(p::Point, l1::Point, l2::Point)
    x, y = l1.coordinates
    dx = l2.coordinates[1] - x
    dy = l2.coordinates[2] - y

    if !iszero(dx) || !iszero(dy)
        t = ((p.coordinates[1] - x) * dx + (p.coordinates[2] - y) * dy) / (dx * dx + dy * dy)

        if t > 1
            x = l2.coordinates[1]
            y = l2.coordinates[2]
        elseif t > 0
            x += dx * t
            y += dy * t
        end
    end

    dx = p.coordinates[1] - x
    dy = p.coordinates[2] - y

    return dx * dx + dy * dy
end


function simplify_radial_distance(points::Vector, tol::Real)
    previous = points[1]
    new_points = [previous]
    point = nothing

    for i in eachindex(points)
        point = points[i]
        if dist²(point, previous) > tol
            push!(new_points, point)
            previous = point
        end
    end

    !isequal(previous, point) && push!(new_points, point)

    return new_points
end

function simplify_DP_step(points::Vector, first::Integer, last::Integer, tol::Real, simplified)
    max_dist = tol
    index = 0

    for i = first+1:last
        dist = segdist²(points[i], points[first], points[last])

        if dist > max_dist
            index = i
            max_dist = dist
        end
    end

    if max_dist > tol
        (index - first > 1) && simplify_DP_step(points, first, index, tol, simplified)
        push!(simplified, points[index])
        (last - index > 1) && simplify_DP_step(points, index, last, tol, simplified)
    end
end


function douglas_peucker(points::Vector, tol::Real)
    last = length(points)

    simplified = [points[1]]
    simplify_DP_step(points, 1, last, tol, simplified)
    push!(simplified, points[last])

    return simplified
end

function DP_simplify(points::Vector, tol::Real, hq::Bool)
    length(points) <= 2 && return points
    sq_tol = tol * tol

    points = hq ? points : simplify_radial_distance(points, sq_tol)
    points = douglas_peucker(points, sq_tol)

    return points
end
