"""
Takes a Point and calculates the ellipse polygon given two semi-axes expressed in variable units and steps for precision.
"""
function ellipse(; center::Point, xAxis::Real, yAxis::Real, steps::Real=64, angle::Real=0., pivot::Union{Point, Nothing}=nothing, units::String="kilometers")
    pivot == nothing && (pivot = center)

    coords = center.coordinates

    if units === "degrees"
        angle = deg2rad(angle)
    else
        xAxis = rhumb_destination(coords, xAxis, 90, units).coordinates
        yAxis = rhumb_destination(coords, yAxis, 0, units).coordinates
        xAxis = xAxis[1] - coords[1]
        yAxis = yAxis[2] - coords[2]
    end

    elCoords = []

    for i in 1:steps
        stepAngle = i * (-360 / steps)
        x = (xAxis * yAxis) / sqrt(yAxis^2 + xAxis^2 * tan(deg2rad(stepAngle))^2)
        y = (xAxis * yAxis) / sqrt(xAxis^2 + yAxis^2 / tan(deg2rad(stepAngle))^2)

        (stepAngle < -90 && stepAngle >= -270) && (x = -x)
        (stepAngle < -180 && stepAngle >= -360) && (y = -y)

        if units === "degrees"
            newx = x * cos(angle) + y * sin(angle)
            newy = y * cos(angle) - x * sin(angle)
            x = newx
            y = newy
        end

        push!(elCoords, [x + coords[1], y + coords[2]])
    end
    push!(elCoords, elCoords[1])

    units === "degrees" && return Polygon([elCoords])

    return transform_rotate(geojson=Polygon([elCoords]), angle=angle, pivot=pivot)

end
