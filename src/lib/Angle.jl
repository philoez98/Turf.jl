"""
    angle_adjacent(start::Position, mid::Position, stop::Position, explementary::Bool=false, mercator::Bool=false)

Find the angle formed by two adjacent segments defined by 3 points. The result will be the (positive clockwise)
angle with origin on the `start-mid` segment, or its explementary angle if required.

# Examples
```julia
julia> using Turf

julia> angle_adjacent(p1, p2, p3, false, false)
38.98929595940595

julia> angle_adjacent(p1, p2, p3, true, false) # explementary angle
321.0107040405941
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
