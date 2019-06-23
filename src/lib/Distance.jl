function distance(from::Position, to::Position, options::String="kilometers")
    δlat = deg2rad((ycoord(to) - ycoord(from)))
    δlon = deg2rad((xcoord(to) - xcoord(from)))

    latFrom = deg2rad(ycoord(from))
    latTo = deg2rad(ycoord(to))

    res = sin(δlat / 2)^2 + sin(δlon / 2)^2 * cos(latFrom) * cos(latTo)

    return radiansToLength(2 * atan(sqrt(res), sqrt(1 - res)), options)
end


function distance(from::Point, to::Point, options::String="kilometers")
    from = from.coordinates
    to = to.coordinates

    δlat = deg2rad((ycoord(to) - ycoord(from)))
    δlon = deg2rad((xcoord(to) - xcoord(from)))

    latFrom = deg2rad(ycoord(from))
    latTo = deg2rad(ycoord(to))

    res = sin(δlat / 2)^2 + sin(δlon / 2)^2 * cos(latFrom) * cos(latTo)

    return radiansToLength(2 * atan(sqrt(res), sqrt(1 - res)), options)
end


function rhumbDistance(from::Position, to::Position, units::String)
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

    res =  Δ * earthRadius

    return convertLength(res, "metres", units)
end

function distanceToSegment(point::Point, first::Point, last::Point, units::String="degrees", method::String="planar")
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


function nearestPoint(target::Point, points::Vector{Point})
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
