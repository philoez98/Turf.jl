"""
    centroid([geojson::T]) where {T <: AbstractGeometry}

Takes one or more features and calculates the centroid using the mean of all vertices.
This lessens the effect of small islands and artifacts when calculating the centroid of a set of polygons.

# Examples
```jldoctest
julia> using Turf

julia> line = LineString([[1, 2], [4, 6], [8, 9.5], [12, 13.4]])
LineString(Array{Float64,1}[[1.0, 2.0], [4.0, 6.0], [8.0, 9.5], [12.0, 13.4]])

julia> centroid(line)
Point([6.25, 7.725])
```
"""
function centroid(geojson::T)::Point where {T<:AbstractGeometry}
    x = 0.
    y = 0.
    len = 0.
    # TODO: support MultiPolygon

    data = geojson.coordinates
    type = geotype(geojson)

    if type === :Point
        x = data[1]
        y = data[2]
        len = 1.
    elseif type === :LineString || type === :MultiPoint
        for i in eachindex(data)
            x += data[i][1]
            y += data[i][2]
            len += 1
        end
    elseif type === :Polygon || type === :MultiLineString
        for i in eachindex(data[1])
            x += data[1][i][1]
            y += data[1][i][2]
            len += 1
        end
    end

    return Point([x / len, y / len])
end

"""
    center(geojson::T) where {T <: AbstractGeometry}

Take a GeoJson Geometry and return the absolute center point.

# Examples
```jldoctest
julia> using Turf

julia> line = LineString([[1, 2], [4, 6], [8, 9.5], [12, 13.4]])
LineString(Array{Float64,1}[[1.0, 2.0], [4.0, 6.0], [8.0, 9.5], [12.0, 13.4]])

julia> center(line)
Point([6.5, 7.7])
```
"""
function center(geojson::T)::Point where {T <: AbstractGeometry}
    box = bbox(geojson)

    x = (box[1] + box[3]) / 2
    y = (box[2] + box[4]) / 2

    return Point([x, y])
end

"""
    masscenter(geojson::T) where {T <: AbstractGeometry}

Take any GeoJson Geometry and return its [center of mass](https://en.wikipedia.org/wiki/Center_of_mass) using this formula:
[Centroid of Polygon](https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon).

# Examples
```jldoctest
julia> using Turf

julia> line = LineString([[1, 2], [4, 6], [8, 9.5], [12, 13.4]])
LineString(Array{Float64,1}[[1.0, 2.0], [4.0, 6.0], [8.0, 9.5], [12.0, 13.4]])

julia> masscenter(line)
Point([6.25, 7.725])
```
"""
function masscenter(geojson::T) where {T <: AbstractGeometry}
    type = geotype(geojson)


    type === :Point && return Point(geojson.coordinates)
    if type === :Polygon
        coords = deepcopy(geojson.coordinates)

        center::Point = centroid(geojson)

        trans = center.coordinates
        sx = 0.
        sy = 0.
        sArea = 0.

        for i in eachindex(coords[1])
            coords[1][i][1] = coords[1][i][1] - trans[1]
            coords[1][i][2] = coords[1][i][2] - trans[2]
        end


        for i in 1:length(coords[1]) - 1
            pi = coords[1][i]
            xi, yi = pi

            pj = coords[1][i + 1]
            xj, yj = pj

            a = xi * yj - xj * yi
            sArea += a

            sx += (xi + xj) * a
            sy += (yi + yj) * a
        end

        sArea == 0 && return center

        area = sArea * 0.5
        factor = 1 / (6 * area)

        return Point([trans[1] + factor * sx, trans[2] + factor * sy])

    else
        # TODO: add convex hull
        return centroid(geojson)

    end
end

"""
    meancenter(geojson::T, weight::Real=1.) where {T <: Union{AbstractGeometry, AbstractFeatureCollection}}

Take a GeoJson Geometry and return the mean center. Can be weighted.

# Examples
```jldoctest
julia> using Turf

julia> line = LineString([[1, 2], [4, 6], [8, 9.5], [12, 13.4]])
LineString(Array{Float64,1}[[1.0, 2.0], [4.0, 6.0], [8.0, 9.5], [12.0, 13.4]])

julia> meancenter(line)
Point([6.25, 7.725])
```
"""
function meancenter(geojson::T, weight::Real=1.) where {T <: Union{AbstractGeometry, AbstractFeatureCollection}}
    # TODO: so ugly! Refactor?
    xs = 0.
    ys = 0.
    ns = 0.

    type = geotype(geojson)
    coors = []

    weight <= 0 && throw(error("Invalid weight. Weight must be > 0."))

    if type === :FeatureCollection

        for feat in geojson.features
            geom = feat.geometry
            coors = geom.coordinates

            if geotype(geom) === :Point
                xs += coors[1] * weight
                ys += coors[2] * weight
                ns += weight

            elseif geotype(geom) === :Polygon || geotype(geom) === :MultiLineString
                for i in eachindex(coors[1])
                    xs += coors[1][i][1] * weight
                    ys += coors[1][i][2] * weight
                    ns += weight
                end
            else
                for i in eachindex(coors)
                    xs += coors[i][1] * weight
                    ys += coors[i][2] * weight
                    ns += weight
                end

            end
        end

    else
        coors = geojson.coordinates

        if type === :Point
            xs += coors[1] * weight
            ys += coors[2] * weight
            ns += weight

        elseif type === :Polygon || type === :MultiLineString
            for i in eachindex(coors[1])
                xs += coors[1][i][1] * weight
                ys += coors[1][i][2] * weight
                ns += weight
            end
        else
            for i in eachindex(coors)
                xs += coors[i][1] * weight
                ys += coors[i][2] * weight
                ns += weight
            end
        end
    end

    return Point([xs / ns, ys / ns])
end

"""
    mediancenter(geojson::T, weight::Real=1., tol::Real=0.001, count::Integer=10) where {T <: AbstractFeatureCollection}

Take a FeatureCollection of points and calculate the median center,
algorithimically. The median center is understood as the point that is
requires the least total travel from all other points.

# Examples
```jldoctest
julia> using Turf

julia> fc = FeatureCollection([Feature(Point([0, 0])),Feature(Point([9, 9])), Feature(Point([9.25, 9.25])),
             Feature(Point([9.5, 9.5])), Feature(Point([9.75, 9.75])), Feature(Point([10, 10]))])
[...]

julia> mediancenter(fc)
Point([9.25425, 9.25425])
```
"""
function mediancenter(geojson::T, weight::Real=1., tol::Real=0.001, count::Integer=10) where {T <: AbstractFeatureCollection}

    mean = meancenter(geojson)

    feat1 = Feature(nothing)

    feat1.properties = Dict(
        "tolerance" => tol,
        "weight" => weight,
        "candidates" => [])

    centroids = FeatureCollection([feat1])

    for feat in geojson.features
        push!(centroids.features, Feature(centroid(feat.geometry)))
    end



    return findMedian(mean.coordinates, [0., 0.], centroids, count)

end

function findMedian(candidate::Position, previous::Position, centroids::FeatureCollection, counter::Integer)

    i = findfirst(x -> x.properties["candidates"] != nothing, centroids.features)
    feat = centroids.features[i]

    tol = feat.properties["tolerance"]
    weight = feat.properties["weight"]
    xSum = 0.
    ySum = 0.
    kSum = 0.
    count = 0

    weight <= 0 && throw(error("Invalid weight. Weight must be > 0."))

    for feats in centroids.features

        feats.geometry == nothing && continue

        count += 1
        dist = weight * distance(feats.geometry, Point(candidate))

        dist === 0 && (dist = 1.)

        k = weight / dist

        xSum += feats.geometry.coordinates[1] * k
        ySum += feats.geometry.coordinates[2] * k

        kSum += k
    end

    count < 1 && throw(error("No features to measure."))

    candidateX = xSum / kSum
    candidateY = ySum / kSum

    (count === 1 || counter === 0 || (abs(candidateX - previous[1])) < tol && (abs(candidateY - previous[2])) < tol) &&
        return Point([candidateX, candidateY])

    push!(feat.properties["candidates"], [candidateX, candidateY])

    return findMedian([candidateX, candidateY], candidate, centroids, counter - 1)
end
