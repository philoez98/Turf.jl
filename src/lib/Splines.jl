import Base.@kwdef

abstract type AbstractSpline end

@kwdef struct Spline{R, P} <: AbstractSpline
    points::Vector{P}
    duration::R = 10000
    sharpness::R = 0.85
    centers::Vector{P} = []
    controls::Vector{Vector{P}} = []
    stepLength::R = 60
    length::R = length(points)
    delay::R = 0
    steps::Vector{R} = []

    function Spline(points::Vector{P}, duration::R, sharpness::R,
        centers::Vector{P},
        controls::Vector{Vector{P}},
        stepLength::R,
        length::R,
        delay::R,
        steps::Vector{R}) where {P, R}

        for point in points
            point.coordinates[3] = length(point) > 2 ? point.coordinates[3] : 0
        end

        for i in eachindex(points)
            p1 = points[i]
            p2 = points[i + 1]

            push!(centers,
                Point([(p1.coordinates[1] + p2.coordinates[1]) / 2,
                (p1.coordinates[2] + p2.coordinates[2]) / 2,
                (p1.coordinates[3] + p2.coordinates[3]) / 2]))
        end
        push!(controls, [points[1], points[1]])

        for j in eachindex(centers)
            c1 = centers[j]
            c2 = centers[j + 1]
            dx = points[j + 1].coordinates[1] - (c1.coordinates[1] + c2.coordinates[1]) / 2
            dy = points[j + 1].coordinates[2] - (c1.coordinates[2] + c2.coordinates[2]) / 2
            dz = points[j + 1].coordinates[3] - (c1.coordinates[3] + c2.coordinates[3]) / 2

            push!(controls, [
                Point([(1. - sharpness) * points[j + 1].coordinates[1] + sharpness * c1.coordinates[1] + dx,
                (1. - sharpness) * points[j + 1].coordinates[2] + sharpness * c1.coordinates[2] + dy,
                (1. - sharpness) * points[j + 1].coordinates[3] + sharpness * c1.coordinates[3] + dz]),
                Point([(1. - sharpness) * points[j + 1].coordinates[1] + sharpness * c2.coordinates[1] + dx,
                (1. - sharpness) * points[j + 1].coordinates[2] + sharpness * c2.coordinates[2] + dy,
                (1. - sharpness) * points[j + 1].coordinates[3] + sharpness * c2.coordinates[3] + dz])
                ])
        end

        push!(controls, [points[length - 1], points[length - 1]])
        steps = cacheSteps(stepLength)

        new{R, P}(points, duration, sharpness, centers, controls, stepLength, length, delay, steps)
    end

end # struct



function cacheSteps(minDist::Real)
    steps = Real[]
    lastStep = pos(0)
    push!(steps, 0)

    for i in 1:10:duration
        step = pos(i)
        dist = sqrt(
            (step.coordinates[1] - lastStep.coordinates[1])^2 +
            (step.coordinates[2] - lastStep.coordinates[2])^2 +
            (step.coordinates[3] - lastStep.coordinates[3])^2)

        if dist > minDist
            push!(steps, i)
            lastStep = step
        end
    end

    return steps

end

function pos(time::Real)
    t = time - delay
    t < 0 && (t = 0)
    t > duration && (t = duration - 1)

    t2 = t / duration
    t2 >= 1 && return points[length - 1]

    n = floor((length(points) - 1) * t2)
    t1 = (length - 1) * t2 - n

    return bezier(t1, points[n], controls[n][1], controls[n + 1][0], points[n + 1])
end

function bezier(; t::Real, p1::Point, c1::Point, c2::Point, p2::Point)
    b = [t^3, (3 * t^2 * (1-t)), (3 * t * (1 - t) * (1 - t)), ((1 - t) * (1 - t) * (1 - t))]

    return Point([
        p2.coordinates[1] * b[1] + c2.coordinates[1] * b[2] + c1.coordinates[1] * b[3] + p1.coordinates[1] * b[4],
        p2.coordinates[2] * b[1] + c2.coordinates[2] * b[2] + c1.coordinates[2] * b[3] + p1.coordinates[2] * b[4],
        p2.coordinates[3] * b[1] + c2.coordinates[3] * b[2] + c1.coordinates[3] * b[3] + p1.coordinates[3] * b[4]])
end
