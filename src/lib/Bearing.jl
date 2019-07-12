"""
    bearing_to_azimuth(bearing::Real)::Real

Convert any bearing angle from the north line direction (positive clockwise)
and return an angle between 0-360 degrees (positive clockwise), 0 being the north line.

# Examples
```julia
julia> using Turf

julia> bearing_to_azimuth(40)
40

julia> bearing_to_azimuth(-105)
255

julia> bearing_to_azimuth(410)
50
```
"""
function bearing_to_azimuth(bearing::Real)::Real
    angle = bearing % 360
    if angle < 0
        angle += 360
    end
    return angle
end


"""
    rhumb_bearing(start::Position, stop::Position, final::Bool=false)

Take two Positions and finds the bearing angle between them along a Rhumb line
i.e. the angle measured in degrees start the north line (0 degrees).

# Examples
```julia
julia> using Turf

julia> start = Position([-75, 45])
2-element Array{Float64,1}:
 -75.0
  45.0

julia> stop = Position([20, 60])
2-element Array{Float64,1}:
 20.0
 60.0

julia> rhumb_bearing(start, stop)
75.28061364784332
```
"""
function rhumb_bearing(start::Position, stop::Position, final::Bool=false)
    bear360 = nothing

    if final === true
        bear360 = calculateRhumbBearing(stop, start)
    else
        bear360 = calculateRhumbBearing(start, stop)
    end

    bear180 = (bear360 > 180) ? - (360 - bear360) : bear360

    return bear180
end

rhumb_bearing(start::Point, stop::Point, final::Bool=false) = rhumb_bearing(start.coordinates, stop.coordinates, final)

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
    bearing(start::Position, stop::Position, final::Bool=false)

Take two Points or Positions and finds the geographic bearing between them,
i.e. the angle measured in degrees from the north line (0 degrees).

# Examples
```julia
julia> using Turf

julia> start = Position([-75, 45])
2-element Array{Float64,1}:
 -75.0
  45.0

julia> stop = Position([20, 60])
2-element Array{Float64,1}:
 20.0
 60.0

julia> bearing(start, stop)
37.75495852601734
```
"""
function bearing(start::Position, stop::Position, final::Bool=false)
    if final === true
        bear = bearing(stop, start, false)
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

bearing(start::Point, stop::Point, final::Bool=false) = bearing(start.coordinates, stop.coordinates, final)
