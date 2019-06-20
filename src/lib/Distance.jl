include("../Utils.jl")
include("../geojson/Geometries.jl")


function distance(from::Geometries.Position, to::Geometries.Position, options::String="kilometers")
    δlat = deg2rad((Geometries.latitude(to) - Geometries.latitude(from)))
    δlon = deg2rad((Geometries.longitude(to) - Geometries.longitude(from)))

    latFrom = deg2rad(Geometries.latitude(from))
    latTo = deg2rad(Geometries.latitude(to))

    res = sin(δlat / 2)^2 + sin(δlon / 2)^2 * cos(latFrom) * cos(latTo)

    return radiansToLength(2 * atan(sqrt(res), sqrt(1 - res)), options)
end


function distance(from::Geometries.Point, to::Geometries.Point, options::String="kilometers")
    from = from.coordinates
    to = to.coordinates

    δlat = deg2rad((Geometries.latitude(to) - Geometries.latitude(from)))
    δlon = deg2rad((Geometries.longitude(to) - Geometries.longitude(from)))

    latFrom = deg2rad(Geometries.latitude(from))
    latTo = deg2rad(Geometries.latitude(to))

    res = sin(δlat / 2)^2 + sin(δlon / 2)^2 * cos(latFrom) * cos(latTo)

    return radiansToLength(2 * atan(sqrt(res), sqrt(1 - res)), options)
end


function rhumbDistance(from::Geometries.Position, to::Geometries.Position, units::String)
    toLon = Geometries.longitude(to)
    toLon += (Geometries.longitude(to) - Geometries.longitude(from) > 180) ? -360 : ((Geometries.longitude(from) - Geometries.longitude(to) > 180)) ? 360 : 0

    ϕ1 = Geometries.latitude(from) * pi / 180
    ϕ2 = Geometries.latitude(to) * pi / 180

    Δϕ = ϕ2 - ϕ1

    Δλ = abs((toLon - Geometries.longitude(from))) * pi / 180
    if Δλ > pi
        Δλ -= 2 * pi
    end

    Δψ = log(tan(ϕ2 / 2 + pi / 4) / tan(ϕ1 / 2 + pi / 4))
    q = abs(Δψ) > 10e-12 ?  Δϕ / Δψ : cos(ϕ1)

    Δ = sqrt(Δϕ * Δϕ + q * q * Δλ * Δλ)

    res =  Δ * Constants.earthRadius

    return convertLength(res, "metres", units)
end

function distanceToSegment(point::Geometries.Point, first::Geometries.Point, last::Geometries.Point, units::String="degrees", method::String="planar")
    v = [last.coordinates[1] - first.coordinates[1], last.coordinates[2], first.coordinates[2]]
    w = [point.coordinates[1] - first.coordinates[1], point.coordinates[2] - first.coordinates[2]]

    c1 = sum(w.*v)
    if c1 <= 0
        return method === "planar" ? rhumbDistance(point.coordinates, first.coordinates, units) : distance(point.coordinates, first.coordinates, "degrees")
    end
    c2 = sum(v.*v)

    if c2 <= c1
        return method === "planar" ? rhumbDistance(point.coordinates, last.coordinates, units) : distance(point.coordinates, last.coordinates, "degrees")
    end

    b2 = c1 / c2

    p = [first.coordinates[1] + (b2 * v[1]), first.coordinates[2] + (b2 * v[2])]

    return method === "planar" ? rhumbDistance(point.coordinates, p.coordinates, units) : distance(point.coordinates, p.coordinates, "degrees")
end


function nearestPoint(target::Geometries.Point, points::Geometries.Points)
    minDistance = Inf
    index = 0

    for (i, point) in enumerate(points)
        dist = distance(target.coordinates, point.coordinates)
        if dist < minDistance
            index = i
            minDistance = dist
        end
    end
    nearest::Geometries.Point = points[index]
    nearest.distance = minDistance

    return nearest
end
