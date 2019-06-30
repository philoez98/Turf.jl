"""Create a Point grid from a bounding box"""
function pointGrid(bbox::Vector{T}, cellSide::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}
    results = []

    west = bbox[1]
    south = bbox[2]
    east = bbox[3]
    north = bbox[4]

    xfrac = cellSide / (distance([west, south], [east, south], units))
    cellWidth = xfrac * (east - west)
    yfrac = cellSide / (distance([west, south], [west, north], units))
    cellHeight = yfrac * (north - south)

    bbWidth = east - west
    bbHeight = north - south

    cols = floor(bbWidth / cellWidth)
    rows = floor(bbHeight / cellHeight)

    Δx  = (bbWidth - cols * cellWidth) / 2
    Δy = (bbHeight - rows * cellHeight) / 2

    curx = west + Δx

    while curx <= east
        cury = south + Δy

        while cury <= north
            pt = Point([curx, cury])
            if mask != nothing
                # if within...
                push!(results, pt)
            else
                push!(results, pt)
            end

            cury += cellHeight
        end

        curx += cellWidth
    end

    return FeatureCollection([Feature(x) for x in results])
end

"""Creates a grid of rectangles from a bounding box."""
function rectangleGrid(bbox::Vector{T}, width::T, height::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}
    results = []

    west = bbox[1]
    south = bbox[2]
    east = bbox[3]
    north = bbox[4]

    xfrac = width / (distance([west, south], [east, south], units))
    cellWidth = xfrac * (east - west)
    yfrac = height / (distance([west, south], [west, north], units))
    cellHeight = yfrac * (north - south)

    bbWidth = east - west
    bbHeight = north - south

    cols = floor(bbWidth / cellWidth)
    rows = floor(bbHeight / cellHeight)

    Δx  = (bbWidth - cols * cellWidth) / 2
    Δy = (bbHeight - rows * cellHeight) / 2

    curx = west + Δx

    for i in 1:cols
        cury = south + Δy

        for j in 1:rows
            poly = Polygon([[
                [curx, cury], [curx, cury + cellHeight], [curx + cellWidth, cury + cellHeight],
                [curx + cellWidth, cury], [curx, cury]]])

            if mask != nothing
                # intersects...
                push!(results, poly)
            else
                push!(results, poly)
            end
            cury += cellHeight
        end
        curx += cellWidth
    end

    return FeatureCollection([Feature(x) for x in results])
end

function hexagon(center::Point, x::T, y::T, cosines::Vector, sines::Vector) where {T <: Real}
    vertices = []
    coords = center.coordinates

    for i in 1:6
        x = coords[1] + x * cosines[i]
        y = coords[2] + y * sines[i]
        push!(vertices, [x, y])
    end
    push!(vertices, vertices[1])

    return Polygon([vertices])
end

function hexTriangles(center::Point, x::T, y::T, cosines::Vector, sines::Vector) where {T <: Real}
    triangles = []
    coords = center.coordinates

    for i in 1:6
        vertices = []
        push!(vertices, coords)
        push!(vertices, [coords[1] + x * cosines[i], coords[2] + y * sines[i]])
        push!(vertices, [coords[1] + x * cosines[(i + 1) % 6], coords[2] + y * sines[(i + 1) % 6]])
        push!(vertices, coords)

        push!(triangles, Polygon([vertices]))
    end

    return triangles
end


"""
Takes a bounding box and the diameter of the cell and returns a FeatureCollection of flat-topped
hexagons or triangles Polygon aligned in an "odd-q" vertical grid as
described in [Hexagonal Grids](http://www.redblobgames.com/grids/hexagons/).
"""
function hexGrid(bbox::Vector{T}, cellSide::T, mask::Union{Polygon, Nothing}=nothing, triangles::Bool=false, units::String="kilometers") where {T <: Real}
    west = bbox[1]
    south = bbox[2]
    east = bbox[3]
    north = bbox[4]

    centerY = (south + north) / 2.
    centerX = (west + east) / 2.

    xfrac = cellSide * 2. / (distance([west, centerY], [east, centerY], units))
    cellWidth = xfrac * (east - west)
    yfrac = cellSide * 2. / (distance([centerX, south], [centerX, north], units))
    cellHeight = yfrac * (north - south)
    radius = cellWidth / 2.

    width = radius * 2.
    height = sqrt(3) / 2. * cellHeight

    bWidth = east - west
    bHeight = north - south

    intx = 3. / 4. * width
    inty = height

    spanx = (bWidth - width) / (width - radius / 2.)
    countx = floor(spanx)

    adjustx = ((countx * intx - radius / 2.) - bWidth) / 2. - radius / 2. + intx / 2.

    county = floor((bHeight - height) / height)
    adjusty = (bHeight - county * height) / 2.

    hasOffsetY = county * height - bHeight > height / 2.

    hasOffsetY && (adjusty -= height / 4.)

    cosines = []
    sines = []

    for i in 1:6
        angle = 2. * pi / 6. * (i - 1)
        push!(cosines, cos(angle))
        push!(sines, sin(angle))
    end

    results = []

    for i = 1:(countx+1)
        for j = 1:(county+1)
            odd = (i-1) % 2 === 1

            #(j === 1 && odd) && continue
            #(j === 1 && hasOffsetY) && continue

            centerx = (i-1) * intx + west - adjustx
            centery = j * inty + south + adjusty

            odd && (centery -= height / 2.)

            if triangles
                t = hexTriangles(Point([centerx, centery]),
                    cellWidth / 2., cellHeight / 2., cosines, sines)

                for triangle in t
                    if mask != nothing
                        # intersects...
                        push!(results, triangle)
                    else
                        push!(results, triangle)
                    end
                end
            else
                hex = hexagon(Point([centerx, centery]), cellWidth / 2., cellHeight / 2., cosines, sines)
                if mask != nothing
                    # intersects...
                    push!(results, hex)
                else
                    push!(results, hex)
                end
            end
        end
    end

    return FeatureCollection([Feature(x) for x in results])
end
