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

    if geotype(geojson) === :Point
        x = data[1]
        y = data[2]
        len = 1
    elseif geotype(geojson) === :LineString
        for i in eachindex(data)
            x += data[i][1]
            y += data[i][2]
            len += 1
        end
    elseif geotype(geojson) === :Polygon
        for i in eachindex(data[1])
            x += data[1][i][1]
            y += data[1][i][2]
            len += 1
        end
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
