include("Constants.jl")

"""
Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function radiansToLength(radians::Number, units::String="kilometers")::Number
    factor = Constants.factors[units]
    if factor == nothing
        throw(error("$(units) is not a valid unit."))
    end
    return radians*factor
end

"""
Convert a distance measurement (assuming a spherical Earth) from a real-world unit to radians.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function lengthToRadians(distance::Number, units::String="kilometers")::Number
    factor =  Constants.factors[units]
    if factor == nothing
        throw(error("$(units) is not a valid unit."))
    end
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
Converts a length to the requested unit.
"""
function convertLength(length::Number, originalUnit::String="kilometers", finalUnit::String="kilometers")::Number
    return length >= 0 ? radiansToLength(lengthToRadians(length, originalUnit), finalUnit) : error("'length' must be a positive number.")
end


"""
Converts an area to the requested unit.
"""
function convertArea(area::Number, originalUnit::String="meters", finalUnit::String="kilometers")::Number
    if area < 0
        throw(error("'area' must be a positive number."))
    end

    startFactor =  Constants.areaFactors[originalUnit]
    endFactor =  Constants.areaFactors[finalUnit]

    return (area / startFactor) * endFactor
end

@inline function allequal(x)
    length(x) < 2 && return true
    e1 = x[1]
    i = 2
    @inbounds for i=2:length(x)
        x[i] == e1 || return false
    end
    return true
end
