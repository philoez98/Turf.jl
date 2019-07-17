"""
    angle_adjacent(start::Position, mid::Position, stop::Position, explementary::Bool=false, mercator::Bool=false)

Find the angle formed by two adjacent segments defined by 3 points. The result will be the (positive clockwise)
angle with origin on the `start-mid` segment, or its explementary angle if required.

# Examples
```jldoctest
julia> p1 = Point([-35, 55])
Point([-35.0, 55.0])

julia> p2 = Point([-34.8, 57.5])
Point([-34.8, 57.5])

julia> p3 = Point([-33.4, 59.1])
Point([-33.4, 59.1])

julia> angle_adjacent(p1, p2, p3, false, false)
202.8279033760424

julia> angle_adjacent(p1, p2, p3, true, false) # explementary angle
157.1720966239576
```
"""
function angle_adjacent(start::Position, mid::Position, stop::Position, explementary::Bool=false, mercator::Bool=false)
    azimuth1 = bearing_to_azimuth((mercator !== true) ? bearing(start, mid, false) : rhumb_bearing(start, mid))
    azimuth2 = bearing_to_azimuth((mercator !== true) ? bearing(stop, mid, false) : rhumb_bearing(stop, mid))

    res = abs(azimuth1 - azimuth2)

    if explementary === true
        return 360 - res
    end

    return res
end

angle_adjacent(start::Point, mid::Point, stop::Point, explementary::Bool=false, mercator::Bool=false) = angle_adjacent(start.coordinates, mid.coordinates, stop.coordinates, explementary, mercator)
