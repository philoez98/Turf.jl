# TODO: implement bbox for featurecollections
# TODO: check correctness
module BBox

include("GeoJSON.jl")
include("Geometries.jl")
include("Features.jl")

using .GeoJSON, .Geometries, .Features

export BBox2D, BBox3D, bbox

const BBox2D = Vector{Float64}(undef, 4)
const BBox3D = Vector{Float64}(undef, 6)


"""
Takes a set of features, calculates the bbox of all input features, and returns a bounding box.
"""
function bbox(geojson::GeoJson)::BBox2D
    result::BBox2D = [Inf, Inf, -Inf, -Inf]

    if geojson.content.type === "FeatureCollection"
        throw(error("Unsupported type"))
    end

    coords = geojson.content.geometry.coordinates

    if typeof(coords) === Vector{Vector{Float64, 1}, 1}

        for (i, el) in coords
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
    end

    return result
end

function bbox(geojson::Geometry)::BBox2D
    result::BBox2D = [Inf, Inf, -Inf, -Inf]

    coords = geojson.coordinates

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

function bbox(geojson::Feature)::BBox2D
    result::BBox2D = [Inf, Inf, -Inf, -Inf]

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

end # module
