using GeoInterface: LineString, Polygon, geotype

function lineSegment(geojson::LineString)
    result::Vector{LineString} = []
    lineSegmentFeature(geojson, result)

    return result
end

function lineSegmentFeature(geojson::Union{LineString, Polygon}, result::Vector{LineString})
    coords = geojson.coordinates

    #geotype(geojson) === :Polygon && (coords = geojson.coordinates)
    #geotype(geojson) === :LineString && (coords = geojson.coordinates)


    segments = createSegments(coords)

    for s in segments
        push!(result, s)
    end

    return result
end

@inline function createSegments(coords::Vector{Position})
    segments = []
    if length(coords) === 2
        push!(segments, Linestring([coords[1], coords[2]]))
        return segments
    end

    for i in 1:length(coords)-1
        line = LineString([coords[i], coords[i + 1]])
        push!(segments, line)
    end
    return segments
end
