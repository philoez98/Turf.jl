module Utils

include("Constants.jl")
include("./geojson/GeoJSON.jl")
include("./geojson/Geometries.jl")

using .GeoJSON

using .Constants: factors, areaFactors, earthRadius
using .Geometries

export radiansToLength, lengthToRadians, lengthToDegrees, bearingToAzimuth, convertLength,
        convertArea, distanceToSegment, rhumbBearing, rhumbDestination, destination,
        bearing, angle


"""
Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function radiansToLength(; radians::Number, units::String="kilometers")::Number
    let factor = factors[units] || throw(error("$(units) is not a valid unit."))
    return radians*factor
end

"""
Convert a distance measurement (assuming a spherical Earth) from a real-world unit to radians.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function lengthToRadians(; distance::Number, units::String="kilometers")::Number
    factor = factors[units] || throw(error("$(units) is not a valid unit."))
    return distance / factor
end

"""
Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
"""
function lengthToDegrees(distance::Number, units::String="kilometers")::Number
    return rad2deg(lengthToRadians(distance, units))
end


"""
Converts any bearing angle from the north line direction (positive clockwise)
and returns an angle between 0-360 degrees (positive clockwise), 0 being the north line
"""
function bearingToAzimuth(bearing::Number)::Number
    angle = bearing % 360
    if angle < 0
        angle += 360
    end
    return angle
end


"""
Converts a length to the requested unit.
"""
function convertLength(length::Float64, originalUnit::String="kilometers", finalUnit::String="kilometers")::Float64
    return length >= 0 ? radiansToLength(lengthToRadians(length, originalUnit), finalUnit) : error("'length' must be a positive number.") end
end


"""
Converts an area to the requested unit.
"""
function convertArea(area::Float64, originalUnit::String="meters", finalUnit::String="kilometers")::Float64
    if area < 0
        throw(error("'area' must be a positive number."))
    end

    startFactor = areaFactors[originalUnit] || error("Invalid units")
    endFactor = areaFactors[finalUnit] || error("Invalid units")

    return (area / startFactor) * endFactor
end


end # module
