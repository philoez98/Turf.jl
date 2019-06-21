using GeoInterface: Position
include("Bearing.jl")

"""
Finds the angle formed by two adjacent segments defined by 3 points. The result will be the (positive clockwise)
angle with origin on the `start-mid` segment, or its explementary angle if required.
"""
function angleAdjacent(start::Position, mid::Position, stop::Position, explementary::Bool, mercator::Bool)
    azimuth1 = bearingToAzimuth((mercator !== true) ? bearing(start, mid, false) : rhumbBearing(start, mid))
    azimuth2 = bearingToAzimuth((mercator !== true) ? bearing(stop, mid, false) : rhumbBearing(stop, mid))

    res = abs(azimuth1 - azimuth2)

    if explementary === true
        return 360 - res
    end

    return res
end
