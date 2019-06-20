include("../geojson/GeoJSON.jl")
include("../geojson/Features.jl")
include("../geojson/Geometries.jl")
include("../geojson/BBox.jl")

"""
    centroid([geojson::GeoJSON.GeoJson])

Takes one or more features and calculates the centroid using the mean of all vertices.
This lessens the effect of small islands and artifacts when calculating the centroid of a set of polygons.
"""
function centroid(geojson::GeoJSON.GeoJson)
    x = 0
    y = 0
    len = 0

    let data end
    if geojson.content.type === "Feature"
        data = geojson.content.geometry.coordinates
    else
        throw(error("FeatureCollections are not supported yet."))
    end

    if typeof(data) === Vector{Vector{Float64, 1}, 1}

        for (i, point) in enumerate(data)
            x += point.coordinates[i][1]
            y += point.coordinates[i][2]
            len += 1
        end
    end

    return Geometries.Point([x / len, y / len])
end

"""Takes a GeoJson and returns the absolute center point of all features."""
function center(geojson::GeoJSON.GeoJson)
    box = bbox(geojson)

    x = (box[1] + box[3]) / 2
    y = (box[2] + box[4]) / 2

    return Geometries.Point([x, y])
end
