
isdefined(Turf, :BBox2D) || const BBox2D = Vector{Float64}(undef, 4)
isdefined(Turf, :BBox3D) || const BBox3D = Vector{Float64}(undef, 6)



"""
    bbox(geojson::T) where {T<:AbstractFeatureCollection}

Take a set of features, calculate the bbox of all input features, and returns a bounding box.
"""
function bbox(geojson::T) where {T<:AbstractFeatureCollection}
    result = [Inf, Inf, -Inf, -Inf]

    for feat in geojson.features
        geom = feat.geometry
        coords = geom.coordinates

        if geotype(geom) === :LineString || geotype(geom) === :MultiPoint
            for i in eachindex(coords)
                result[1] > coords[i][1] && (result[1] = coords[i][1])

                result[2] > coords[i][2] && (result[2] = coords[i][2])

                result[3] < coords[i][1] && (result[3] = coords[i][1])

                result[4] < coords[i][2] && (result[4] = coords[i][2])
            end
        elseif geotype(geom) === :Polygon || geotype(geom) === :MultiLineString
            for j in eachindex(coords)
                for i in eachindex(coords[j])
                    result[1] > coords[j][i][1] && (result[1] = coords[j][i][1])

                    result[2] > coords[j][i][2] && (result[2] = coords[j][i][2])

                    result[3] < coords[j][i][1] && (result[3] = coords[j][i][1])

                    result[4] < coords[j][i][2] && (result[4] = coords[j][i][2])
                end
            end
        elseif geotype(geom) === :Point
                result[1] > coords[1] && (result[1] = coords[1])

                result[2] > coords[2] && (result[2] = coords[2])

                result[3] < coords[1] && (result[3] = coords[1])

                result[4] < coords[2] && (result[4] = coords[2])
        elseif geotype(geom) === :MultiPolygon
            for i in eachindex(coords)
                for j in eachindex(coords[i])
                    for k in eachindex(coords[i][j])
                        result[1] > coords[i][j][k][1] && (result[1] = coords[i][j][k][1])

                        result[2] > coords[i][j][k][2] && (result[2] = coords[i][j][k][2])

                        result[3] < coords[i][j][k][1] && (result[3] = coords[i][j][k][1])

                        result[4] < coords[i][j][k][2] && (result[4] = coords[i][j][k][2])
                    end
                end
            end
        end
    end
    return result
end

"""
    bbox(geojson::T) where {T <: AbstractGeometry}

Take a GeoJSON Geometry and calculate its bounding box.
"""
function bbox(geojson::T) where {T <: AbstractGeometry}
    result = [Inf, Inf, -Inf, -Inf]

    coords = geojson.coordinates

    if geotype(geojson) === :LineString || geotype(geojson) === :MultiPoint
        for i in eachindex(coords)
            result[1] > coords[i][1] && (result[1] = coords[i][1])

            result[2] > coords[i][2] && (result[2] = coords[i][2])

            result[3] < coords[i][1] && (result[3] = coords[i][1])

            result[4] < coords[i][2] && (result[4] = coords[i][2])
        end
    elseif geotype(geojson) === :Polygon || geotype(geojson) === :MultiLineString
        for j in eachindex(coords)
            for i in eachindex(coords[j])
                result[1] > coords[j][i][1] && (result[1] = coords[j][i][1])

                result[2] > coords[j][i][2] && (result[2] = coords[j][i][2])

                result[3] < coords[j][i][1] && (result[3] = coords[j][i][1])

                result[4] < coords[j][i][2] && (result[4] = coords[j][i][2])
            end
        end
    elseif geotype(geojson) === :Point
            result[1] > coords[1] && (result[1] = coords[1])

            result[2] > coords[2] && (result[2] = coords[2])

            result[3] < coords[1] && (result[3] = coords[1])

            result[4] < coords[2] && (result[4] = coords[2])
    elseif geotype(geojson) === :MultiPolygon
        for i in eachindex(coords)
            for j in eachindex(coords[i])
                for k in eachindex(coords[i][j])
                    result[1] > coords[i][j][k][1] && (result[1] = coords[i][j][k][1])

                    result[2] > coords[i][j][k][2] && (result[2] = coords[i][j][k][2])

                    result[3] < coords[i][j][k][1] && (result[3] = coords[i][j][k][1])

                    result[4] < coords[i][j][k][2] && (result[4] = coords[i][j][k][2])
                end
            end
        end
    end

    return result
end

"""
     bbox(geojson::T) where {T<: AbstractFeature}

Take a Feature and return a bounding box around its geometry.
"""
bbox(geojson::T) where {T<: AbstractFeature} = bbox(geojson.geometry)


"""
    bbox_polygon(bbox::Vector{T}) where {T <: Real}

Take a bbox and return an equivalent Polygon.
"""
function bbox_polygon(bbox::Vector{T})::Polygon where {T <: Real}
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


function clip_to_bbox(geom::AbstractGeometry, box::Vector{T}) where {T <: Real}
    type = geotype(geom)
    coords = geom.coordinates

    if isequal(type, :LineString) || isequal(type, :MultiLineString)
        lines = []

        type === :Linestring && (coords = [geom.coordinates])
        println(coords)
        for line in coords
            println(lineclip([Point(line)], box))
            push!(lines, lineclip([Point(line)], box))
        end

        length(lines) === 1 && return LineString(lines[1])
        return MultiLineString(lines)
    elseif isequal(type, :Polygon)
        return Polygon(polyclip(coords, box))
    elseif isequal(type, :MultiPolygon)
        return MultiPolygon(map(x -> return polyclip(x, box), coords))
    else
        throw(error("Unsupported Geometry."))
    end
end

function polyclip(rings, box::Vector{T}) where {T <: Real}
    out = []
    for ring in rings
        clipped = polygonclip([Point(x) for x in ring], box)

        if length(clipped) > 0
            if (clipped[1][1] !== clipped[end][1] || clipped[1][2] !== clipped[end][2])
                push!(clipped, clipped[1])
            end
            length(clipped) >= 4 && push!(out, clipped)
        end
    end

    return out
end
