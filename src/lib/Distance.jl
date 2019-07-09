
"""
    distance(from::Position, to::Position, units::String="kilometers")

Calculate the distance between two Points or Positions in degrees, radians, miles, or kilometers.
It use the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
"""
function distance(from::Position, to::Position, units::String="kilometers")
    δlat = deg2rad((ycoord(to) - ycoord(from)))
    δlon = deg2rad((xcoord(to) - xcoord(from)))

    latFrom = deg2rad(ycoord(from))
    latTo = deg2rad(ycoord(to))

    res = sin(δlat / 2)^2 + sin(δlon / 2)^2 * cos(latFrom) * cos(latTo)

    return radians_to_length(2 * atan(sqrt(res), sqrt(1 - res)), units)
end


distance(from::Point, to::Point, units::String="kilometers") = distance(from.coordinates, to.coordinates, units)

"""
    rhumb_distance(from::Position, to::Position, units::String="kilometers")

Calculate the distance along a rhumb line between two Points or Positions in degrees, radians,
miles, or kilometers.
"""
function rhumb_distance(from::Position, to::Position, units::String="kilometers")
    toLon = xcoord(to)
    toLon += (xcoord(to) - xcoord(from) > 180) ? -360 : ((xcoord(from) - xcoord(to) > 180)) ? 360 : 0

    ϕ1 = ycoord(from) * pi / 180
    ϕ2 = ycoord(to) * pi / 180

    Δϕ = ϕ2 - ϕ1

    Δλ = abs((toLon - xcoord(from))) * pi / 180
    if Δλ > pi
        Δλ -= 2 * pi
    end

    Δψ = log(tan(ϕ2 / 2 + pi / 4) / tan(ϕ1 / 2 + pi / 4))
    q = abs(Δψ) > 10e-12 ?  Δϕ / Δψ : cos(ϕ1)

    Δ = sqrt(Δϕ * Δϕ + q * q * Δλ * Δλ)

    res =  Δ * earth_radius

    return convert_length(res, "metres", units)
end

rhumb_distance(from::Point, to::Point, units::String="kilometers") = rhumb_distance(from.coordinates, to.coordinates, units)


function distance_to_segment(point::Point, first::Point, last::Point, units::String="degrees", method::String="planar")
    v = [last.coordinates[1] - first.coordinates[1], last.coordinates[2], first.coordinates[2]]
    w = [point.coordinates[1] - first.coordinates[1], point.coordinates[2] - first.coordinates[2]]

    c1 = sum(w.*v)
    if c1 <= 0
        return method === "planar" ? rhumb_distance(point.coordinates, first.coordinates, units) : distance(point.coordinates, first.coordinates, "degrees")
    end
    c2 = sum(v.*v)

    if c2 <= c1
        return method === "planar" ? rhumb_distance(point.coordinates, last.coordinates, units) : distance(point.coordinates, last.coordinates, "degrees")
    end

    b2 = c1 / c2

    p = [first.coordinates[1] + (b2 * v[1]), first.coordinates[2] + (b2 * v[2])]

    return method === "planar" ? rhumb_distance(point.coordinates, p.coordinates, units) : distance(point.coordinates, p.coordinates, "degrees")
end

"""
    nearestpoint(target::Point, points::Vector{Point})

Take a reference Point and an array of Points and return the point from the
array closest to the reference. This calculation is geodesic.
"""
function nearestpoint(target::Point, points::Vector{Point})
    minDistance = Inf
    index = 0

    for (i, point) in enumerate(points)
        dist = distance(target.coordinates, point.coordinates)
        if dist < minDistance
            index = i
            minDistance = dist
        end
    end
    nearest::Point = points[index]

    return nearest
end

"""
    distance_weight(; geojson::T, treshold::Real=10000, p::Real=2, binary::Bool=false, alpha::Real=-1., standardization::Bool=false) where {T <: AbstractFeatureCollection}

Calculate the Minkowski p-norm distance between two Points.
"""
function distance_weight(; geojson::T, treshold::Real=10000, p::Real=2, binary::Bool=false, alpha::Real=-1., standardization::Bool=false) where {T <: AbstractFeatureCollection}
    features = []

    for feat in geojson.features
        push!(features, centroid(feat.geometry))
    end

    weights = []
    for i in 1:length(geojson.features)
        weights[i] = []
    end

    for i in 1:length(geojson.features)
        for j in i:length(geojson.features)

            i === j && (weights[i][j] = 0.)

            dist = pnorm_distance(geojson.features[i], geojson.features[j], p)
            weights[i][j] = dist
            weights[j][i] = dist
        end
    end

    for i in 1:length(geojson.features)
        for j in 1:length(geojson.features)
            dist = weights[i][j]
            dist === 0 && continue

            if binary
                if dist <= treshold
                    weights[i][j] = 1.
                else
                    weights[i][j] = 0.
                end
            else
                if dist <= treshold
                    weights[i][j] = dist^alpha
                else
                    weights[i][j] = 0.
                end

            end
        end
    end

    if standardization
        for i in 1:length(geojson.features)
            rows = reduce(+, weights[i])

            for j in eachindex(geojson.features)
                weights[i][j] = weights[i][j] / rows
            end
        end
    end

    return weights
end

"""
    pnorm_distance(point1::Point, point2::Point, p::Real=2)::Real

Calculate the Minkowski p-norm distance between two Points.
"""
function pnorm_distance(point1::Point, point2::Point, p::Real=2)::Real
    coords1 = point1.coordinates
    coords2 = point2.coordinates

    Δx = coords1[1] - coords2[1]
    Δy = coords1[2] - coords2[2]

    p === 1 && return abs(Δx) + abs(Δy)

    return (Δx^p + Δy^p)^(1/p)

end

"""
    midpoint(first::Point, second::Point)

Take two Points and returns a point midway between them.
The midpoint is calculated geodesically, meaning the curvature of the earth is taken into account.
"""
function midpoint(first::Point, second::Point)
    dist = distance(first, second)
    bear = bearing(first.coordinates, second.coordinates)

    return destination(first.coordinates, dist / 2, bear)
end
