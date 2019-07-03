"""
Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function radians_to_length(radians::Number, units::String="kilometers")::Number
    factor = factors[units]
    if factor == nothing
        throw(error("$(units) is not a valid unit."))
    end
    return radians*factor
end

"""
Convert a distance measurement (assuming a spherical Earth) from a real-world unit to radians.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function length_to_radians(distance::Number, units::String="kilometers")::Number
    factor =  factors[units]
    if factor == nothing
        throw(error("$(units) is not a valid unit."))
    end
    return distance / factor
end

"""
Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
"""
function length_to_degrees(distance::Number, units::String="kilometers")::Number
    return rad2deg(length_to_radians(distance, units))
end

"""
Converts a length to the requested unit.
"""
function convert_length(length::Number, originalUnit::String="kilometers", finalUnit::String="kilometers")::Number
    return length >= 0 ? radians_to_length(length_to_radians(length, originalUnit), finalUnit) : error("'length' must be a positive number.")
end


"""
Converts an area to the requested unit.
"""
function convert_area(area::Number, originalUnit::String="meters", finalUnit::String="kilometers")::Number
    if area < 0
        throw(error("'area' must be a positive number."))
    end

    startFactor =  area_factors[originalUnit]
    endFactor =  area_factors[finalUnit]

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

"""Convert 900913 x/y values to lon/lat."""
function to_WGS84(pos::Point)
    a = 6378137.0

    return [rad2deg(pos.coordinates[1]) / a,
        rad2deg((pi * 0.5) - atan(exp(-pos.coordinates[2] / a)))]
end

"""Convert lon/lat values to 900913 x/y."""
function to_mercator(pos::Point)
    a = 6378137.0
    extent = 20037508.342789244

    coords = pos.coordinates

    adjusted = (abs(coords[1]) <= 180) ? coords[1] : (coords[1] - (sign(coords[1] * 360)))
    xy = [deg2rad(a * adjusted), a * log(atan((pi * 0.25) + (deg2rad(0.5 * coords[2]))))]

    xy[1] > extent && (xy[1] = extent)
    xy[1] < -extent && (xy[1] = -extent)
    xy[2] > extent && (xy[2] = extent)
    xy[2] < -extent && (xy[2] = -extent)

    return xy
end

"""Takes one or more features and returns their area in square meters."""
function area(geojson::T) where {T <: Union{AbstractFeatureCollection, AbstractGeometry}}
    if geotype(geojson) === :FeatureCollection
        results = []
        for feat in geojson.features
            push!(results, calculateArea(feat.geometry))
        end

        return results
    else
        return calculateArea(geojson)
    end
end

function calculateArea(geom::T) where {T <: AbstractGeometry}
    total = 0
    type = geotype(geom)

    (type === :Point || type === :LineString || type === :MultiLineString) && return 0

    type === :Polygon && return polygonArea(geom.coordinates)

    if type === :MultiPolygon
        for i in 1:length(geom.coordinates)
            total += polygonArea(geom.coordinates[i])
        end

        return total
    end

    return 0
end

function polygonArea(coords)
    total = 0

    if length(coords) > 0
        total += abs(ringArea(coords[1]))

        for i in 2:length(coords)
            total -= abs(ringArea(coords[i]))
        end
    end

    return total
end


function ringArea(coords)
    total = 0
    lowIndex = 1
    midIndex = 1
    upIndex = 1

    if length(coords) > 2
        for i in 1:length(coords) - 1
            if i === length(coords) - 2
                lowIndex = length(coords) - 2
                midIndex = length(coords) - 1
                upIndex = 1
            elseif i === length(coords) - 1
                lowIndex = length(coords) - 1
                midIndex = 1
                upIndex = 2
            else
                lowIndex = i
                midIndex = i + 1
                upIndex = i + 2
            end

            p1 = coords[lowIndex]
            p2 = coords[midIndex]
            p3 = coords[upIndex]

            total += (p3[1] * (pi / 180) - (p1[1] * (pi / 180))) * sin(p2[2] * (pi / 180))
        end

        total = total * 6378137 * 6378137 / 2
    end

    return total
end
