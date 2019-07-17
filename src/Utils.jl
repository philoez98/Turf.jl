"""
    radians_to_length(radians::Real, units::String="kilometers")::Real

Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function radians_to_length(radians::Real, units::String="kilometers")::Real
    factor = factors[units]
    if factor == nothing
        throw(error("$(units) is not a valid unit."))
    end
    return radians*factor
end

"""
    length_to_radians(distance::Real, units::String="kilometers")::Real

Convert a distance measurement (assuming a spherical Earth) from a real-world unit to radians.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
"""
function length_to_radians(distance::Real, units::String="kilometers")::Real
    factor =  factors[units]
    if factor == nothing
        throw(error("$(units) is not a valid unit."))
    end
    return distance / factor
end

"""
    length_to_degrees(distance::Real, units::String="kilometers")::Real

Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees.
Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
"""
function length_to_degrees(distance::Real, units::String="kilometers")::Real
    return rad2deg(length_to_radians(distance, units))
end

"""
    convert_length(length::Real, originalUnit::String="kilometers", finalUnit::String="kilometers")::Real

Convert a length to the requested unit.
"""
function convert_length(length::Real, originalUnit::String="kilometers", finalUnit::String="kilometers")::Real
    return length >= 0 ? radians_to_length(length_to_radians(length, originalUnit), finalUnit) : error("'length' must be a positive number.")
end


"""
    convert_area(area::Real, originalUnit::String="meters", finalUnit::String="kilometers")::Real

Convert an area to the requested unit.
"""
function convert_area(area::Real, originalUnit::String="meters", finalUnit::String="kilometers")::Real
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

"""
    to_WGS84(pos::Point)

Convert 900913 x/y values to lon/lat.
"""
function to_WGS84(pos::Point)
    a = 6378137.0

    return [rad2deg(pos.coordinates[1]) / a,
        rad2deg((pi * 0.5) - atan(exp(-pos.coordinates[2] / a)))]
end

"""
    to_mercator(pos::Point)

Convert lon/lat values to 900913 x/y.
"""
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

"""
    area(geojson::T) where {T <: Union{AbstractFeatureCollection, AbstractGeometry}}
Take one or more features and return their area in square meters.
"""
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

"""
    clean(geojson::AbstractGeometry, mutate::Bool=false)

Remove redundant coordinates from any GeoJSON Geometry.
"""
function clean(geojson::AbstractGeometry, mutate::Bool=false)

    type = geotype(geojson)
    coords = []

    isequal(type, :Point) && return geojson

    !mutate && (geojson = deepcopy(geojson))

    if isequal(type, :LineString) || isequal(type, :MultiPoint)
        coords = cleanline(geojson)
    elseif isequal(type, :MultiLineString) || isequal(type, :Polygon)
        foreach(x -> push!(coords, cleanline(x)), geojson.coordinates[1])
    elseif isequal(type, :MultiPolygon)
        points = []

        foreach(x -> push!(points, cleanline(x)), geojson.coordinates)
        push!(coords, points)
    else
        throw(error("$(type) is not supported."))

    end

    mutate && (geojson.coordinates = coords); return geojson
    return Feature(geojson)
end

clean!(geojson::AbstractGeometry) = clean(geojson, true)

function cleanline(line)
    points = line.coordinates

    (length(points) == 2 && point[1] != point[2]) && return points

    new_points = []
    index = length(points) - 1
    np_length = length(new_points)

    push!(new_points, points[1])

    for i in 2:index
        prev = new_points[length(new_points)]

        (points[i][1] == prev[1] && points[i][2] == prev[2]) && continue

        push!(new_points, points[i])
        np_length = length(new_points)

        if np_length > 2
            if point_on_line(Point(new_points[np_length - 2]), LineString([new_points[np_length], new_points[np_length- 1]]))
                splice!(new_points, length(new_points)-1)
            end
        end
    end

    push!(new_points, points[length(points)])

    np_length = length(new_points)

    (points[1] == points[end] && np_length < 4) && throw(error("invalid polygon."))
    if point_on_line(Point(new_points[np_length - 2]), LineString([new_points[np_length], new_points[np_length- 1]]))
        splice!(new_points, length(new_points)-1)
    end

    return new_points
end
