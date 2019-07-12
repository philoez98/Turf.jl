"""
    destination(origin::Position, distance::Real, bearing::Real, units::String="kilometers")

Take a Point or a Position and calculate the location of a destination point given a distance in
degrees, radians, miles, or kilometers; and bearing in degrees.
The destination is calculated using the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.

# Examples
```julia
julia> using Turf

julia> point = Point([-75, 38])
Point([-75.0, 38.0])

julia> destination(point, 100, 0)
Point([-75.0, 38.8993])

julia> destination(point, 100, 45)
Point([-74.1859, 38.6331])

julia> destination(point, 50, 0, "miles")
Point([-75.0, 38.7237])
```
"""
function destination(origin::Position, distance::Real, bearing::Real, units::String="kilometers")
    lon1 = deg2rad(origin[1])
    lat1 = deg2rad(origin[2])
    bearingRad = deg2rad(bearing)
    radians = length_to_radians(distance, units)

    lat2 = asin(sin(lat1) * cos(radians) + cos(lat1) * sin(radians) * cos(bearingRad))
    lon2 = lon1 + atan(sin(bearingRad) * sin(radians) * cos(lat1), cos(radians) - sin(lat1) * sin(lat2))

    return Point([rad2deg(lon2), rad2deg(lat2)])
end

destination(origin::Point, distance::Real, bearing::Real, units::String="kilometers") = destination(origin.coordinates, distance, bearing, units)


"""
    rhumb_destination(origin::Position, distance::Real, bearing::Real, units::String="kilometers")

Take a Point or a Position and return the destination Point having travelled the given distance along a Rhumb line from the
origin Point with the (varant) given bearing.

# Examples
```julia
julia> using Turf

julia> point = Point([-75, 38])
Point([-75.0, 38.0])

julia> rhumb_destination(point, 100, 0)
Point([-75.0, 38.8993])

julia> rhumb_destination(point, 100, 45)
Point([-74.1895, 38.6359])

julia> rhumb_destination(point, 50, 0, "miles")
Point([-75.0, 38.7237])
```
"""
function rhumb_destination(origin::Position, distance::Real, bearing::Real, units::String="kilometers")
    negative::Bool = distance < 0
    distanceinMeters = convert_length(abs(distance), units, "meters")
    if negative
        distanceinMeters = - abs(distanceinMeters)
    end

    dest = calculateRhumbDestination(origin, distanceinMeters, bearing)

    dest[1] += (dest[1] - origin[1] > 180) ? -360 : (origin[1] - dest[1] > 180) ? 360 : 0

    return Point(dest)
end

rhumb_destination(origin::Point, distance::Real, bearing::Real, units::String="kilometers") = rhumb_destination(origin.coordinates, distance, bearing, units)

function calculateRhumbDestination(origin::Position, distance::Real, bearing::Real)
    Δ = distance / earth_radius
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
