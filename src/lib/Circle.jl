"""
    circle(; center::Union{Point, Position}, radius::Real=5., steps::Integer=64, units::String="kilometers")

Take a Point or a Position and calculate the circle polygon given a radius in degrees, radians, miles, or kilometers; and steps for precision.

# Examples
```julia
julia> using Turf

julia> point = Point([35, 45])
Point([35.0, 45.0])

julia> circle(center=point, steps=5)
Polygon(Array{Array{Float64,1},1}[[[34.9395, 45.0139], [34.9626, 44.9636], [35.0374, 44.9636], [35.0605, 45.0139], [35.0, 45.045], [34.9395, 45.0139]]])

julia> circle(center=point, radius=2.5, steps=5, units="degrees")
Polygon(Array{Array{Float64,1},1}[[[31.5893, 45.7231], [32.99, 42.9571], [37.01, 42.9571], [38.4107, 45.7231], [35.0, 47.5029], [31.5893, 45.7231]]])
```
"""
function circle(; center::Union{Point, Position}, radius::Real=5., steps::Integer=64, units::String="kilometers")
    coords = []
    position::Position = []

    if geotype(center) === :Point
        position = Position(center.coordinates)
    else
        position = center
    end

    for i in 1:steps
        push!(coords, destination(position, Float64(radius), i * (-360 / steps), units).coordinates)
    end
    push!(coords, coords[1])

    return Polygon([coords])
end

"""
    sector(center::Point, radius::Real, bearing1::Real, bearing2::Real, steps::Integer=64, units::String="kilometers")

Creates a circular sector of a circle of given radius and center Point,
between (clockwise) bearing1 and bearing2; 0 bearing is North of center point, positive clockwise.

# Examples
```julia
julia> using Turf

julia> point = Point([35, 45])
Point([35.0, 45.0])

julia> sector(point, 5, 0, 0.5, 5)
Polygon(Array{Array{Float64,1},1}[[[35.0, 45.0], [35.0, 45.045], [35.0006, 45.045], [35.0, 45.0]]])

julia> sector(point, 5, 0, 270, 5)
Polygon(Array{Array{Float64,1},1}[[[35.0, 45.0], [35.0, 45.045], [35.0605, 45.0139], [35.0374, 44.9636], [34.9626, 44.9636], [34.9364, 45.0], [35.0, 45.0]]])
```
"""
function sector(center::Point, radius::Real, bearing1::Real, bearing2::Real, steps::Integer=64, units::String="kilometers")
    complem(bearing1) === complem(bearing2) && return circle(center=center, radius=radius, steps=steps, units=units)

    coords = center.coordinates
    arc = linearc(center, radius, bearing1, bearing2, steps, units)
    sliceCoords = [[coords]]

    for i in eachindex(arc.coordinates)
        push!(sliceCoords[1], arc.coordinates[i])
    end
    push!(sliceCoords[1], coords)

    return Polygon(sliceCoords)
end

function complem(a::Real)
    b = a % 360
    b < 0 && (b += 360)

    return b
end
