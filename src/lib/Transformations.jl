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

function explode(geojson::T)::FeatureCollection where {T <: Union{AbstractFeatureCollection, AbstractGeometry}}
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

    return FeatureCollection([Feature(x) for x in points])

end
