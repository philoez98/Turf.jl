module Splines

include("../geojson/Geometries.jl")

import Base.@kwdef

abstract type AbstractSpline end

@kwdef struct Spline <: AbstractSpline
    points::Geometries.Points
    duration::Real = 10000
    sharpness::Real = 0.85
    centers::Geometries.Points = []
    controls::Vector{Geometries.Points} = []
    stepLength::Real = 60
    length::Real = length(points)
    delay::Real = 0
    steps::Vector{Real} = []

    function Spline(points::Geometries.Points, duration::Real, sharpness::Real,
                    centers::Geometries.Points,
                    controls::Vector{Geometries.Points},
                    stepLength::Real,
                    length::Real,
                    delay::Real,
                    steps::Vector{Real})

        for i in 1:length(points)
            p1 = points[i]
            p2 = points[i + 1]

            push!(centers,
            Geometries.Point([(p1.coordinates[1] + p2.coordinates[1]) / 2,
            (p1.coordinates[2] + p2.coordinates[2]) / 2,
            (p1.coordinates[3] + p2.coordinates[3]) / 2]))
        end
        push!(controls, [points[1], points[1]])

        for j in 1:length(centers)
            c1 = centers[j]
            c2 = centers[j + 1]
            dx = points[j + 1].coordinates[1] - (c1.coordinates[1] + c2.coordinates[1]) / 2
            dy = points[j + 1].coordinates[2] - (c1.coordinates[2] + c2.coordinates[2]) / 2
            dz = points[j + 1].coordinates[3] - (c1.coordinates[3] + c2.coordinates[3]) / 2

            push!(controls, [
            Geometries.Point([(1. - sharpness) * points[j + 1].coordinates[1] + sharpness * c1.coordinates[1] + dx,
            (1. - sharpness) * points[j + 1].coordinates[2] + sharpness * c1.coordinates[2] + dy,
            (1. - sharpness) * points[j + 1].coordinates[3] + sharpness * c1.coordinates[3] + dz]),
            Geometries.Point([(1. - sharpness) * points[j + 1].coordinates[1] + sharpness * c2.coordinates[1] + dx,
            (1. - sharpness) * points[j + 1].coordinates[2] + sharpness * c2.coordinates[2] + dy,
            (1. - sharpness) * points[j + 1].coordinates[3] + sharpness * c2.coordinates[3] + dz])
            ])
        end

        push!(controls, [points[length - 1], points[length - 1]])
        steps = cacheSteps(stepLength)
    end
end # struct



function cacheSteps(minDist::Real)

end

end # module
