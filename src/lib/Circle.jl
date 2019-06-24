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
