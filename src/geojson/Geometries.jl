### Based on the GeoInterface.jl implementation
module Geometries

export AbstractGeometry, AbstractPoint, AbstractMultiPoint, AbstractLineString,
        AbstractMultiLineString, AbstractPolygon, AbstractMultiPolygon,
        Position, latitude, longitude, elevation, coordinates, Point, MultiPoint,
        MultiPolygon, MultiLineString, Polygon, LineString, Geometry


const Position = Vector{Float64}
const allowed_types = ["Point", "MultiPoint", "LineString",
      "MultiLineString", "Polygon", "MultiPolygon"]

abstract type AbstractGeometry end
abstract type AbstractPoint <: AbstractGeometry end
abstract type AbstractMultiPoint <: AbstractGeometry end
abstract type AbstractLineString <: AbstractGeometry end
abstract type AbstractMultiLineString <: AbstractGeometry end
abstract type AbstractPolygon <: AbstractGeometry end
abstract type AbstractMultiPolygon <: AbstractGeometry end


latitude(p::Position) = p[2]
longitude(p::Position) = p[1]
elevation(p::Position) = length(p) === 3 ? p[3] : 0


mutable struct Point <: AbstractPoint
    coordinates::Position

    Point(coordinates) = length(coordinates) > 3 || length(coordinates) === 0 ? error("Out of Range") : new(coordinates)
end

coordinates(obj::Vector{T}) where {T <: AbstractPoint} = Position[map(coordinates, obj)...]
coordinates(obj::Vector{T}) where {T <: AbstractLineString} = Vector{Position}[map(coordinates, obj)...]
coordinates(obj::Vector{Vector{T}}) where {T <: AbstractPoint} = Vector{Position}[map(coordinates, obj)...]
coordinates(obj::Vector{Position}) = obj
coordinates(obj::Vector{Vector{Position}}) = obj
coordinates(obj::Vector{Vector{Vector{Position}}}) = obj
coordinates(obj::Vector{Vector{Vector{T}}}) where {T <: AbstractPoint} = Vector{Vector{Position}}[map(coordinates, obj)...]
coordinates(obj::Vector{Vector{T}}) where {T <: AbstractLineString} = Vector{Vector{Position}}[map(coordinates, obj)...]
coordinates(obj::Vector{T}) where {T <: AbstractPolygon} = Vector{Vector{Position}}[map(coordinates, obj)...]

Point(x::Float64,y::Float64) = Point([x,y])
Point(x::Float64,y::Float64,z::Float64) = Point([x,y,z])
Point(point::AbstractPoint) = Point(coordinates(point))

mutable struct MultiPoint <: AbstractMultiPoint
    coordinates::Vector{Position}
end

MultiPoint(point::Position) = MultiPoint(Position[point])
MultiPoint(point::AbstractPoint) = MultiPoint(Position[coordinates(point)])

MultiPoint(points::Vector{T}) where {T <: AbstractPoint} = MultiPoint(coordinates(points))
MultiPoint(points::AbstractMultiPoint) = MultiPoint(coordinates(points))
MultiPoint(line::AbstractLineString) = MultiPoint(coordinates(line))

mutable struct LineString <: AbstractLineString
    coordinates::Vector{Position}
end

LineString(points::Vector{T}) where {T <: AbstractPoint} = LineString(coordinates(points))
LineString(points::AbstractMultiPoint) = LineString(coordinates(points))
LineString(line::AbstractLineString) = LineString(coordinates(line))

mutable struct MultiLineString <: AbstractMultiLineString
    coordinates::Vector{Vector{Position}}
end

MultiLineString(line::Vector{Position}) = MultiLineString(Vector{Position}[line])
MultiLineString(line::Vector{T}) where {T <: AbstractPoint} = MultiLineString(Vector{Position}[coordinates(line)])
MultiLineString(line::AbstractLineString) = MultiLineString(Vector{Position}[coordinates(line)])

MultiLineString(lines::Vector{Vector{T}}) where {T <: AbstractPoint} = MultiLineString(coordinates(lines))
MultiLineString(lines::Vector{T}) where {T <: AbstractLineString} = MultiLineString(Vector{Position}[map(coordinates,lines)])
MultiLineString(lines::AbstractMultiLineString) = MultiLineString(coordinates(lines))
MultiLineString(poly::AbstractPolygon) = MultiLineString(coordinates(poly))

mutable struct Polygon <: AbstractPolygon
    coordinates::Vector{Vector{Position}}
end

Polygon(line::Vector{Position}) = Polygon(Vector{Position}[line])
Polygon(line::Vector{T}) where {T <: AbstractPoint} = Polygon(Vector{Position}[coordinates(line)])
Polygon(line::AbstractLineString) = Polygon(Vector{Position}[coordinates(line)])

Polygon(lines::Vector{Vector{T}}) where {T <: AbstractPoint} = Polygon(coordinates(lines))
Polygon(lines::Vector{T}) where {T <: AbstractLineString} = Polygon(coordinates(lines))
Polygon(lines::AbstractMultiLineString) = Polygon(coordinates(lines))
Polygon(poly::AbstractPolygon) = Polygon(coordinates(poly))

mutable struct MultiPolygon <: AbstractPolygon
    coordinates::Vector{Vector{Vector{Position}}}
end

MultiPolygon(line::Vector{Position}) = MultiPolygon(Vector{Vector{Position}}[Vector{Position}[line]])
MultiPolygon(line::Vector{T}) where {T <: AbstractPoint} = MultiPolygon(Vector{Vector{Position}}[Vector{Position}[coordinates(line)]])
MultiPolygon(line::AbstractLineString) = MultiPolygon(Vector{Vector{Position}}[Vector{Position}[coordinates(line)]])

MultiPolygon(poly::Vector{Vector{T}}) where {T <: AbstractPoint} = MultiPolygon(Vector{Vector{Position}}[coordinates(poly)])
MultiPolygon(poly::Vector{T}) where {T <: AbstractLineString} = MultiPolygon(Vector{Vector{Position}}[coordinates(poly)])
MultiPolygon(poly::AbstractMultiLineString) = MultiPolygon(Vector{Vector{Position}}[coordinates(poly)])
MultiPolygon(poly::AbstractPolygon) = MultiPolygon(Vector{Vector{Position}}[coordinates(poly)])

MultiPolygon(polys::Vector{Vector{Vector{T}}}) where {T <: AbstractPoint} = MultiPolygon(coordinates(polys))
MultiPolygon(polys::Vector{Vector{T}}) where {T <: AbstractLineString} = MultiPolygon(coordinates(polys))
MultiPolygon(polys::Vector{T}) where {T <: AbstractPolygon} = MultiPolygon(coordinates(polys))
MultiPolygon(polys::AbstractMultiPolygon) = MultiPolygon(coordinates(polys))

mutable struct Geometry <: AbstractGeometry
    type::String
    coordinates::Union{Point, LineString, MultiPoint, MultiLineString, Polygon, MultiPolygon}

    Geometry(type::String, coordinates::Union) = type in allowed_types ? new(type, coordinates) : throw(error("Invalid 'type' found."))
end

end # module
