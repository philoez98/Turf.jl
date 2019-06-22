using GeoInterface: LineString, Polygon, geotype

function lineSegment(geojson::LineString)
    result::Vector{LineString} = []
    lineSegmentFeature(geojson, result)

    return result
end

function lineSegmentFeature(geojson::Union{LineString, Polygon}, result::Vector{LineString})
    coords::Vector{Position} = []
    geom = geojson.geometry

    geotype(geom) === :Polygon && (coords = geom.coordinates)
    geotype(geom) === :LineString && (coords = [geom.coordinates])


    segments = createSegments(coords)

    for s in segments
        push!(result, s)
    end

    return result
end

@inline function createSegments(coords::Vector{Position})
    segments = []
    for i in eachindex(coords) - 1
        line = LineString([coords[i], coords[i + 1]])
        push!(segments, line)
    end
end
