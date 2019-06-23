"""
    centroid([geojson::GeoJSON.GeoJson])

Takes one or more features and calculates the centroid using the mean of all vertices.
This lessens the effect of small islands and artifacts when calculating the centroid of a set of polygons.
"""
function centroid(geojson::T) where {T<:AbstractGeometry}
    x = 0
    y = 0
    len = 0
    # TODO: check dims
    data = geojson.coordinates

    for (i, point) in enumerate(data)
        x += point.coordinates[i][1]
        y += point.coordinates[i][2]
        len += 1
    end

    return Point([x / len, y / len])
end

"""Takes a GeoJson and returns the absolute center point of all features."""
function center(geojson::T) where {T <: AbstractGeometry}
    box = bbox(geojson)

    x = (box[1] + box[3]) / 2
    y = (box[2] + box[4]) / 2

    return Point([x, y])
end
