using GeoInterface: Position, Point
include("../Utils.jl")

"""
Takes a Point and calculates the location of a destination point given a distance in
degrees, radians, miles, or kilometers; and bearing in degrees.
This uses the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
"""
function destination(origin::Position, distance::Real, bearing::Real, units::String="kilometers")
    lon1 = deg2rad(origin[1])
    lat1 = deg2rad(origin[2])
    bearingRad = deg2rad(bearing)
    radians = lengthToRadians(distance, units)

    lat2 = asin(sin(lat1) * cos(radians) + cos(lat1) * sin(radians) * cos(bearingRad))
    lon2 = lon1 + atan(sin(bearingRad) * sin(radians) * cos(lat1), cos(radians) - sin(lat1) * sin(lat2))

    return Point([rad2deg(lon2), rad2deg(lat2)])
end


"""
Returns the destination Point having travelled the given distance along a Rhumb line from the
origin Point with the (varant) given bearing.
"""
function rhumbDestination(origin::Position, distance::Real, bearing::Real, units::String="kilometers")
    negative::Bool = distance < 0
    distanceinMeters = convertLength(abs(distance), units, "meters")
    if negative
        distanceinMeters = - abs(distanceinMeters)
    end

    dest = calculateRhumbDestination(origin, distanceinMeters, bearing)

    dest[1] += (dest[1] - origin[1] > 180) ? -360 : (origin[1] - dest[1] > 180) ? 360 : 0

    return Point(dest)
end

function calculateRhumbDestination(origin::Position, distance::Real, bearing::Real)
    Δ = distance / Constants.earthRadius
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
