# TODO: implement bbox for featurecollections
# TODO: check correctness

isdefined(Turf, :BBox2D) || const BBox2D = Vector{Float64}(undef, 4)
isdefined(Turf, :BBox3D) || const BBox3D = Vector{Float64}(undef, 6)



""" Takes a set of features, calculates the bbox of all input features, and returns a bounding box."""
function bbox(geojson::T) where {T<:AbstractFeatureCollection}
    results::Vector{Vector{Float64}} = []
    for feat in geojson.features
        geom = feat.geometry
        push!(results, bbox(geom))
    end

    return results
end

function bbox(geojson::T) where {T <: AbstractGeometry}
    result = [Inf, Inf, -Inf, -Inf]

    coords = geojson.coordinates

    if geotype(geojson) === :LineString
        for i in eachindex(coords)
            result[1] > coords[i][1] && (result[1] = coords[i][1])

            result[2] > coords[i][2] && (result[2] = coords[i][2])

            result[3] < coords[i][1] && (result[3] = coords[i][1])

            result[4] < coords[i][2] && (result[4] = coords[i][2])
        end
    elseif geotype(geojson) === :Polygon
        for i in eachindex(coords[1])
            result[1] > coords[1][i][1] && (result[1] = coords[1][i][1])

            result[2] > coords[1][i][2] && (result[2] = coords[1][i][2])

            result[3] < coords[1][i][1] && (result[3] = coords[1][i][1])

            result[4] < coords[1][i][2] && (result[4] = coords[1][i][2])
        end
    else
            result[1] > coords[1] && (result[1] = coords[1])

            result[2] > coords[2] && (result[2] = coords[2])

            result[3] < coords[1] && (result[3] = coords[1])

            result[4] < coords[2] && (result[4] = coords[2])
    end

    return result
end

function bbox(geojson::T) where {T<: AbstractFeature}
    return bbox(geojson.geometry)
end

"""Takes a bbox and returns an equivalent Polygon."""
function bbox_polygon(bbox::Vector{T}) where {T <: Real}
    west = bbox[1]
    south = bbox[2]
    east = bbox[3]
    north = bbox[4]

    length(bbox) === 6 && throw(error("BBoxes with 6 positions are not supported."))

    lowLeft = [west, south]
    topLeft = [west, north]
    topRight = [east, north]
    lowRight = [east, south]

    return Polygon([[lowLeft, lowRight, topRight, topLeft, lowLeft]])

end
