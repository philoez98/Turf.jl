# Turf

This document tracks the currently implemented features of the Turf.js library.

## Measurement
- [ ] turf-along
- [x] turf-area
- [x] turf-bbox
- [x] turf-bbox-polygon
- [x] turf-bearing
- [x] turf-center
- [x] turf-center-mean
- [x] turf-center-median
- [x] turf-center-of-mass
- [x] turf-centroid
- [x] turf-destination
- [x] turf-distance
- [x] turf-distance-weight
- [ ] turf-envelope
- [ ] turf-length
- [x] turf-midpoint
- [ ] turf-point-on-feature
- [x] turf-polygon-tangents
- [ ] turf-point-to-line-distance
- [x] turf-rhumb-bearing
- [x] turf-rhumb-destination
- [x] turf-rhumb-distance
- [x] turf-square
- [ ] turf-great-circle

## Coordinate Mutation
- [ ] turf-clean-coords
- [x] turf-flip
- [ ] turf-rewind
- [ ] turf-truncate

## Transformation
- [ ] turf-bbox-clip
- [x] turf-bezier-spline
- [ ] turf-buffer
- [x] turf-circle
- [x] turf-ellipse
- [ ] turf-clone
- [x] turf-concave
- [ ] turf-convex
- [ ] turf-difference
- [ ] turf-directional-mean
- [ ] turf-dissolve
- [ ] turf-intersect
- [ ] turf-line-offset
- [ ] turf-simplify
- [ ] turf-tesselate
- [x] turf-transform-rotate
- [x] turf-transform-translate
- [x] turf-transform-scale
- [ ] turf-union
- [ ] turf-voronoi

## Feature Conversion
- [x] turf-combine
- [x] turf-explode
- [ ] turf-flatten
- [ ] turf-line-to-polygon
- [ ] turf-polygonize
- [x] turf-polygon-to-line

## Misc
- [ ] turf-kinks
- [x] turf-line-arc
- [ ] turf-line-chunk
- [ ] turf-line-intersect
- [ ] turf-line-overlap
- [x] turf-line-segment
- [x] turf-line-slice
- [x] turf-line-slice-along
- [ ] turf-line-split
- [ ] turf-mask
- [ ] turf-nearest-point-on-line
- [ ] turf-nearest-point-to-line
- [x] turf-sector
- [ ] turf-points-within-polygon
- [ ] turf-polygon-smooth
- [ ] turf-moran-index
- [ ] turf-invariant
- [ ] turf-shortest-path
- [ ] turf-unkink-polygon

## Helper
The helper functions are all part of the GeoJson module; this class contains unit conversion and additional functionality to help other calculations.

## Random
- [ ] turf-random-interpolate
- [ ] turf-random-isobands
- [ ] turf-random-isolines
- [x] turf-random-planepoint
- [ ] turf-random
- [ ] turf-random-tin

## Data
- [ ] turf-sample

## Joins
- [ ] turf-inside
- [x] turf-tag
- [x] turf-points-within-polygon

## Grids
- [x] turf-hex-grid
- [x] turf-point-grid
- [x] turf-rectangle-grid
- [x] turf-square-grid
- [x] turf-triangle-grid

## Classification
- [x] turf-nearest-point
- [x] turf-nearest-neighbour-analysis
- [ ] turf-quadrat-analysis

## Aggregation
- [ ] turf-collect
- [ ] turf-clusters-dbscan
- [ ] turf-clusters-kmeans
- [ ] turf-clusters

## Meta
The meta functions are not implemented because of the *GeoInterface.jl* module. It already contains all the specific methods to handle geographic data.

## Booleans
- [x] turf-boolean-clockwise
- [x] turf-boolean-contains
- [x] turf-boolean-crosses
- [x] turf-boolean-disjoint
- [x] turf-boolean-equal
- [ ] turf-boolean-intersects
- [ ] turf-boolean-overlap
- [x] turf-boolean-parallel
- [x] turf-boolean-point-in-polygon
- [x] turf-boolean-point-on-line
- [ ] turf-boolean-touches
- [ ] turf-boolean-valid
- [x] turf-boolean-within

## Unit Conversion
- [x] bearingToAzimuth
- [x] convertArea
- [x] convertLength
- [x] degreesToRadians
- [x] lengthToRadians
- [x] lengthToDegrees
- [x] radiansToLength
- [x] radiansToDegrees
- [x] toMercator
- [x] toWgs84

## Others
- [x] turf-projection
