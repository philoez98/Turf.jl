using GeoInterface: AbstractGeometry, Feature,
AbstractFeatureCollection, AbstractFeature, geotype, Point
include("Bearing.jl")
include("BBox.jl")
include("Distance.jl")
include("Centering.jl")


const origin_options = ["sw", "se", "nw", "ne", "center", "centroid"]

"""
Rotates any geojson Feature or Geometry of a specified angle, around its `centroid` or a given `pivot` point;
all rotations follow the right-hand rule.
"""
function transformRotate(geojson::T, angle::Float64) where {T <: AbstractGeometry}
    if angle === 0
        return geojson
    end

    pivot = centroid(geojson)

    for point in geojson.coordinates
        initAngle = rhumbBearing(pivot, point.coordinates)
        finalAngle = initAngle + angle
        dist = rhumbDistance(pivot, point.coordinates)
        newCoords = rhumbDestination(pivot, dist, finalAngle).coordinates
        point.coordinates[1] = newCoords[1]
        point.coordinates[2] = newCoords[2]
    end

    return geojson
end



"""
Scale a GeoJSON from a given point by a factor of scaling (ex: factor=2 would make the GeoJSON 200% larger).
If a FeatureCollection is provided, the origin point will be calculated based on each individual Feature.
"""
function transformScale(geojson::T, factor::Float64, origin::String="centroid") where {T <: AbstractFeatureCollection}

    for (i, feat) in enumerate(geojson.features)
        geojson.features[i] = scale(feat, factor, origin)
    end

    return geojson
end

transformScale(geojson::T, factor::Float64, origin::String="centroid") where {T <: AbstractFeature} = scale(geojson, factor, origin)


"""
Scales a Feature.
"""
function scale(feature::Feature, factor::Float64, origin::String="centroid")

    if !(origin in origin_options)
        throw(error("'$(origin)' is not a valid option. Allowed values are: $(origin_options)"))
    end

    (geotype(feature.geometry) === :Point || factor === 1) && return feature

    coords = feature.geometry.coordinates

    start = rhumbDistance(center, coords)
    bearing = rhumbBearing(center, coord)
    distance = start * factor

    newCoords = rhumbDestination(center, distance, bearing).coordinates

    coords[1] = newCoords[1]
    coords[2] = newCoords[2]

    if length(coords) === 3
        coords[3] *= factor
    end

    return feature

end

function getOrigin(geojson::Feature, origin::String)
    if origin === nothing || origin === ""
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
        return center(geojson)
    elseif origin === "centroid"
        return centroid(geojson)
    else
        throw(error("Invalid origin."))
    end
end
