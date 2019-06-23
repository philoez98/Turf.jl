# TODO: implement bbox for featurecollections
# TODO: check correctness

using GeoInterface: AbstractFeature, AbstractGeometry, AbstractFeatureCollection, geotype

const BBox2D = Vector{Float64}(undef, 4)
const BBox3D = Vector{Float64}(undef, 6)


""" Takes a set of features, calculates the bbox of all input features, and returns a bounding box."""
function bbox(geojson::T) where {T<:AbstractFeatureCollection}
    result = [Inf, Inf, -Inf, -Inf]

    coords = []

    for feat in geojson.features
        push!(coords, feat.geometry.coordinates)
    end


    for (i, el) in enumerate(coords)
        if result[1] > el[i][1]
            result[1] = el[i][1]
        end

        if result[2] > el[i][2]
            result[2] = el[i][2]
        end

        if result[3] < el[i][1]
            result[3] = el[i][1]
        end

        if result[4] < el[i][2]
            result[4] = el[i][2]
        end
    end

    return result
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
        for el in coords
            result[1] > el[1] && (result[1] = el[1])

            result[2] > el[2] && (result[2] = el[2])

            result[3] < el[1] && (result[3] = el[1])

            result[4] < el[2] && (result[4] = el[2])
        end
    end

    return result
end

function bbox(geojson::T) where {T<: AbstractFeature}
    result = [Inf, Inf, -Inf, -Inf]

    coords = geojson.geometry.coordinates

    for el in coords
        if result[1] > el[1]
            result[1] = el[1]
        end

        if result[2] > el[2]
            result[2] = el[2]
        end

        if result[3] < el[1]
            result[3] = el[1]
        end

        if result[4] < el[2]
            result[4] = el[2]
        end
    end

    return result
end
