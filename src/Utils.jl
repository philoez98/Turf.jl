module Utils

include("Constants.jl")
include("./geojson/geom/Geometries.jl")

using .Constants: factors, areaFactors, earthRadius
using .Geometries

export radiansToLength, lengthToRadians, lengthToDegrees, bearingToAzimuth, convertLength,
        convertArea, distanceToSegment, rhumbBearing, rhumbDestination, destination, distance,
        bearing, angle


"""
Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function radiansToLength(; radians::Number, units::String="kilometers")::Number
    let factor = factors[units] || throw(error("$(units) is not a valid unit."))
    return radians*factor
end

"""
Convert a distance measurement (assuming a spherical Earth) from a real-world unit to radians.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function lengthToRadians(; distance::Number, units::String="kilometers")::Number
    factor = factors[units] || throw(error("$(units) is not a valid unit."))
    return distance / factor
end

"""
Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
"""
function lengthToDegrees(distance::Number, units::String="kilometers")::Number
    return rad2deg(lengthToRadians(distance, units))
end


"""
Converts any bearing angle from the north line direction (positive clockwise)
and returns an angle between 0-360 degrees (positive clockwise), 0 being the north line
"""
function bearingToAzimuth(bearing::Number)::Number
    angle = bearing % 360
    if angle < 0
        angle += 360
    end
    return angle
end


"""
Converts a length to the requested unit.
"""
function convertLength(length::Float64, originalUnit::String="kilometers", finalUnit::String="kilometers")::Float64
    return length >= 0 ? radiansToLength(lengthToRadians(length, originalUnit), finalUnit) : error("'length' must be a positive number.") end
end


"""
Converts an area to the requested unit.
"""
function convertArea(area::Float64, originalUnit::String="meters", finalUnit::String="kilometers")::Float64
    if area < 0
        throw(error("'area' must be a positive number."))
    end

    startFactor = areaFactors[originalUnit] || error("Invalid units")
    endFactor = areaFactors[finalUnit] || error("Invalid units")

    return (area / startFactor) * endFactor
end

function distance(from::Position, to::Position, options::String)
    δlat = deg2rad((latitude(to) - latitude(from)))
    δlon = deg2rad((longitude(to) - longitude(from)))

    latFrom = deg2rad(latitude(from))
    latTo = deg2rad(latitude(to))

    res = sin(δlat / 2)^2 + sin(δlon / 2)^2 * cos(latFrom) * cos(latTo)

    return radiansToLength(2 * atan(sqrt(res), sqrt(1 - res)), options)
end

function rhumbDistance(from::Position, to::Position, units::String)
    toLon = longitude(to)
    toLon += (longitude(to) - longitude(from) > 180) ? -360 : ((longitude(from) - longitude(to) > 180)) ? 360 : 0

    ϕ1 = latitude(from) * pi / 180
    ϕ2 = latitude(to) * pi / 180

    Δϕ = ϕ2 - ϕ1

    Δλ = abs((toLon - longitude(from))) * pi / 180
    if Δλ > pi
        Δλ -= 2 * pi
    end

    Δψ = log(tan(ϕ2 / 2 + pi / 4) / tan(ϕ1 / 2 + pi / 4))
    q = abs(Δψ) > 10e-12 ?  Δϕ / Δψ : cos(ϕ1)

    Δ = sqrt(Δϕ * Δϕ + q * q * Δλ * Δλ)

    return Δ * radius
end


"""
Takes two Positions and finds the bearing angle between them along a Rhumb line
i.e. the angle measured in degrees start the north line (0 degrees)
"""
function rhumbBearing(start::Position, stop::Position, final::Bool)
    bear360 = nothing

    if final === true
        bear360 = calculateRhumbBearing(stop, start)
    else
        bear360 = calculateRhumbBearing(start, stop)
    end

    bear180 = (bear360 > 180) ? - (360 - bear360) : bear360

    return bear180

end

function calculateRhumbBearing(a::Position, b::Position)
    ϕ1 = deg2rad(a[2])
    ϕ2 = deg2rad(b[2])

    Δλ = deg2rad((b[1] - a[1]))

    if Δλ > pi Δλ -= 2 * pi end
    if Δλ < -pi Δλ += 2* pi end

    Δψ = log(tan(ϕ2 / 2 + pi / 4) / tan(ϕ1 / 2 + pi / 4))

    θ = atan(Δλ, Δψ)

    return (rad2deg(θ) + 360) % 360
end

"""
Takes a Point and calculates the location of a destination point given a distance in
degrees, radians, miles, or kilometers; and bearing in degrees.
This uses the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
"""
function destination(origin::Position, distance::Float64, bearing::Float64, units::String)
    lon1 = deg2rad(origin[1])
    lat1 = deg2rad(origin[2])
    bearingRad = deg2rad(bearing)
    radians = lengthToRadians(distance, units)

    lat2 = asin(sin(lat1) * cos(radians) + cos(lat1) * sin(radians) * cos(bearingRad))
    lon2 = lon1 + atan(sin(bearingRad) * sin(radians) * cos(lat1) * cos(radians) - sin(lat1) * sin(lat2))

    return Point([deg2rad(lon2), deg2rad(lat2)])
end


"""
Returns the destination Point having travelled the given distance along a Rhumb line from the
origin Point with the (varant) given bearing.
"""
function rhumbDestination(origin::Position, distance::Float64, bearing::Float64, units::String)
    negative::Bool = distance < 0
    distanceinMeters = convertLength(abs(distance), units, "meters")
    if negative
        distanceinMeters = - abs(distanceinMeters)
    end

    dest = calculateRhumbDestination(origin, distanceinMeters, bearing)

    dest[1] += (dest[1] - origin[1] > 180) ? -360 : (origin[1] - dest[1] > 180) ? 360 : 0

    return Point(dest)
end

function calculateRhumbDestination(origin::Position, distance::Float64, bearing::Float64)
    Δ = distance / earthRadius
    λ1 = origin[1] * pi /180
    ϕ1 = deg2rad(origin[2])
    θ = deg2rad(bearing)

    Δϕ = Δ * cos(θ)
    ϕ2 = ϕ1 + Δϕ

    if abs(ϕ2 > pi / 2)
        ϕ2 = ϕ2 > 0 ? pi - ϕ2 : -pi - ϕ2
    end

    Δψ = log(tan(ϕ2 / 2 + pi / 4) / tan(ϕ1 / 2 + pi / 4))
    q = abs(Δψ) > 10e-12 ? Δϕ / Δψ : cos(ϕ1)

    Δλ = Δ * sin(θ) / q
    λ2 = λ1 + Δλ

    return [((λ2 * 180 / pi) + 540) % 360 - 180, ϕ2 * 180 / pi]
end


"""
Takes two points and finds the geographic bearing between them,
i.e. the angle measured in degrees from the north line (0 degrees)
"""
function bearing(start::Position, stop::Position, final::Bool)
    if final === true
        bear = bearing(stop, start)
        return (bear + 180) % 360
    end

    lon1 = deg2rad(start[1])
    lon2 = deg2rad(stop[1])
    lat1 = deg2rad(start[2])
    lat2 = deg2rad(stop[2])

    a = sin(lon2 - lon1) * cos(lat2)
    b = cos(lat1) * sin(lat2) - sin(lat1)  * cos(lat2) * cos(lon2 - lon1)

    return rad2deg(atan(a, b))

end

"""
Finds the angle formed by two adjacent segments defined by 3 points. The result will be the (positive clockwise)
angle with origin on the `start-mid` segment, or its explementary angle if required.
"""
function angle(start::Position, mid::Position, stop::Position, explementary::Bool, mearcator::Bool)
    azimuth1 = bearingToAzimuth((mercator !== true) ? bearing(start, mid) : rhumbBearing(start, mid))
    azimuth2 = bearingToAzimuth((mercator !== true) ? bearing(stop, mid) : rhumbBearing(stop, mid))

    res = abs(azimuth1 - azimuth2)

    if explementary === true
        return 360 - res
    end

    return res
end


"""
Takes one or more features and calculates the centroid using the mean of all vertices.
This lessens the effect of small islands and artifacts when calculating the centroid of a set of polygons.
"""
function centroid(geojson::GeoJson)
    x = 0, y = 0, len = 0

    for point in geojson.geometry
        x += point.coordinates[1]
        y += point.coordinates[2]
        len++
    end

    return Point([x / len, y / len])
end


"""
Rotates any geojson Feature or Geometry of a specified angle, around its `centroid` or a given `pivot` point;
all rotations follow the right-hand rule.
"""
function transformRotate(geojson::GeoJson, angle::Float64)
    if angle === 0 return geojson end

    pivot = centroid(geojson)

    for point in geojson.geometry
        initAngle = rhumbBearing(pivot, point.coordinates)
        finalAngle = initAngle + angle
        dist = rhumbDistance(pivot, point.coordinates)
        newCoords = rhumbDestination(pivot, dist, finalAngle).coordinates
        point.coordinates[1] = newCoords[1]
        point.coordinates[2] = newCoords[2]
    end

    return geojson
end


function nearestPoint(target::Point, points::Vector{Point})
    minDistance = Inf
    index = 0

    for (i, point) in enumerate(points)
        dist = distance(target, point)
        if dist < minDistance
            index = i
            minDistance = dist
        end
    end
    nearest::Point = points[i]
    nearest.distance = minDistance

    return nearest
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

    return method === "planar" ? rhumbDistance(point.coordinates, p.coordinates, units) : distance(point, p, "degrees")
end

end # module
