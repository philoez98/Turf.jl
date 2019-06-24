module Turf

using GeoInterface: Point, Position, LineString, Feature, AbstractGeometry, AbstractFeature,
    AbstractFeatureCollection, xcoord, ycoord, Polygon, geotype, MultiPolygon, MultiLineString, FeatureCollection

export radiansToLength, lengthToRadians, lengthToDegrees, convertLength, convertArea,
    angleAdjacent, distance, distanceToSegment, rhumbDistance, rhumbBearing, rhumbDestination,
    destination, bearing, bearingToAzimuth, bbox, center, centroid, transformRotate, transformScale,
    bezier, concave, clockwise, nearestPoint, earthRadius, areaFactors, unitsFactor, parallel,
    Spline, lineSegment, pointOnLine, pointInPolygon, medianCenter, massCenter, meanCenter


include("Constants.jl")
include("Utils.jl")
include("lib/Angle.jl")
include("lib/Distance.jl")
include("lib/Destination.jl")
include("lib/Bearing.jl")
include("lib/BBox.jl")
include("lib/Lines.jl")
include("lib/Centering.jl")
include("lib/Transformations.jl")
include("lib/Splines.jl")
include("lib/BezierSpline.jl")
include("lib/Booleans.jl")




end
