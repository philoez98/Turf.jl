include("Lines.jl")
"""
Takes a Point or a Position and calculates the circle polygon given a radius in degrees, radians, miles, or kilometers; and steps for precision.
"""
function circle(; center::Union{Point, Position}, radius::Real=5, steps::Integer=64, units::String="kilometers")
    coords = []
    position::Position = []

    if geotype(center) === :Point
        position = Position(center.coordinates)
    else
        position = center
    end

    for i in 1:steps
        push!(coords, destination(position, Float64(radius), i * (-360 / steps)).coordinates)
    end
    push!(coords, coords[1])

    return Polygon([coords])
end

"""
Creates a circular sector of a circle of given radius and center Point,
between (clockwise) bearing1 and bearing2; 0 bearing is North of center point, positive clockwise.
"""
function sector(center::Point, radius::Real, bearing1::Real, bearing2::Real, steps::Real=64., units::String="kilometers")
    to360(bearing1) === to360(bearing2) && return circle(center=center, radius=radius, steps=steps, units=units)

    coords = center.coordinates
    arc = linearc(center, radius, bearing1, bearing2, steps, units)
    sliceCoords = [[coords]]

    for i in eachindex(arc.coordinates)
        push!(sliceCoords[1], arc.coordinates[i])
    end
    push!(sliceCoords[1], coords)

    return Polygon(sliceCoords)
end
