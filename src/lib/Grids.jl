"""
    point_grid(bbox::Vector{T}, cellSide::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}

Create a Point grid from a bounding box

# Examples
```jldoctest
julia> bbox = [-1, 2, 1, 3]
4-element Array{Int64,1}:
 -1
  2
  1
  3

julia> point_grid(bbox, 100)
FeatureCollection{Feature}(Feature[Feature(Point([-0.899869, 2.05034]), Dict{String,Any}()), Feature(Point([-0.899869, 2.94966]),
 Dict{String,Any}()), Feature(Point([0.0, 2.05034]), Dict{String,Any}()), Feature(Point([0.0, 2.94966]), Dict{String,Any}()),
 Feature(Point([0.899869, 2.05034]), Dict{String,Any}()), Feature(Point([0.899869, 2.94966]), Dict{String,Any}())], nothing, nothing
```
"""
function point_grid(bbox::Vector{T}, cellSide::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}
    results = []

    west = bbox[1]
    south = bbox[2]
    east = bbox[3]
    north = bbox[4]

    xfrac = cellSide / (distance(float([west, south]), float([east, south]), units))
    cellWidth = xfrac * (east - west)
    yfrac = cellSide / (distance(float([west, south]), float([west, north]), units))
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

"""
    rectangle_grid(bbox::Vector{T}, width::T, height::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}

Create a grid of rectangles from a bounding box.

# Examples
```jldoctest
julia> bbox = [-1, 2, 1, 3]
4-element Array{Int64,1}:
 -1
  2
  1
  3

julia> rectangle_grid(bbox, 150, 100)
FeatureCollection{Feature}(Feature[Feature(Polygon(Array{Array{Float64,1},1}[[[-0.674901, 2.05034], [-0.674901, 2.94966], [0.674901, 2.94966], [0.674901, 2.05034], [-0.674901, 2.05034]]]), Dict{String,Any}())], nothing, nothing)
```
"""
function rectangle_grid(bbox::Vector{T}, width::T, height::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}
    results = []

    west = bbox[1]
    south = bbox[2]
    east = bbox[3]
    north = bbox[4]

    xfrac = width / (distance(float([west, south]), float([east, south]), units))
    cellWidth = xfrac * (east - west)
    yfrac = height / (distance(float([west, south]), float([west, north]), units))
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
     hexgrid(bbox::Vector{T}, cellSide::T, mask::Union{Polygon, Nothing}=nothing, triangles::Bool=false, units::String="kilometers") where {T <: Real}

Take a bounding box and the diameter of the cell and return a FeatureCollection of flat-topped
hexagons or triangles Polygon aligned in an "odd-q" vertical grid as
described in [Hexagonal Grids](http://www.redblobgames.com/grids/hexagons/).
"""
function hexgrid(bbox::Vector{T}, cellSide::T, mask::Union{Polygon, Nothing}=nothing, triangles::Bool=false, units::String="kilometers") where {T <: Real}
    west = bbox[1]
    south = bbox[2]
    east = bbox[3]
    north = bbox[4]

    centerY = (south + north) / 2.
    centerX = (west + east) / 2.

    xfrac = cellSide * 2. / (distance(float([west, centerY]), float([east, centerY]), units))
    cellWidth = xfrac * (east - west)
    yfrac = cellSide * 2. / (distance(float([centerX, south]), float([centerX, north]), units))
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

            (j === 1 && odd) && continue
            (j === 1 && hasOffsetY) && continue

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

"""
    square_grid(bbox::Vector{T}, cell_side::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}

Create a square grid from a bounding box.

# Examples
```jldoctest
julia> bbox = [-1, 2, 1, 3]
4-element Array{Int64,1}:
 -1
  2
  1
  3

  julia> square_grid(bbox, 100)
  FeatureCollection{Feature}(Feature[Feature(Polygon(Array{Array{Float64,1},1}[[[-0.899869, 2.05034], [-0.899869, 2.94966], [0.0, 2.94966], [0.0, 2.05034], [-0.899869, 2.05034]]]), Dict{String,Any}()),
  Feature(Polygon(Array{Array{Float64,1},1}[[[0.0, 2.05034], [0.0, 2.94966], [0.899869, 2.94966], [0.899869, 2.05034], [0.0, 2.05034]]]), Dict{String,Any}())], nothing, nothing)
```
"""
function square_grid(bbox::Vector{T}, cell_side::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}
    return rectangle_grid(bbox, cell_side, cell_side, mask, units)
end

"""
    triangle_grid(bbox::Vector{T}, cell_side::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}

Take a bounding box and a cell depth and returns a set of triangular Polygons in a grid.

# Examples
```jldoctest
julia> bbox = [-1, 2, 1, 3]
4-element Array{Int64,1}:
 -1
  2
  1
  3

  julia> triangle_grid(bbox, 200)
  FeatureCollection{Feature}(Feature[Feature(Polygon(Array{Array{Float64,1},1}[[[-1.0, 2.0], [-1.0, 3.79864], [0.799737, 2.0], [-1.0, 2.0]]]), Dict{String,Any}()),
  Feature(Polygon(Array{Array{Float64,1},1}[[[-1.0, 3.79864], [0.799737, 3.79864], [0.799737, 2.0], [-1.0, 3.79864]]]), Dict{String,Any}()),
  Feature(Polygon(Array{Array{Float64,1},1}[[[0.799737, 2.0], [0.799737, 3.79864], [2.59947, 3.79864], [0.799737, 2.0]]]), Dict{String,Any}()),
  Feature(Polygon(Array{Array{Float64,1},1}[[[0.799737, 2.0], [2.59947, 3.79864], [2.59947, 2.0], [0.799737, 2.0]]]), Dict{String,Any}())], nothing, nothing)
```
"""
function triangle_grid(bbox::Vector{T}, cell_side::T, mask::Union{Polygon, Nothing}=nothing, units::String="kilometers") where {T <: Real}
    results = []

    x_frac = cell_side / (distance(float([bbox[1], bbox[2]]), float([bbox[3], bbox[2]]), units))
    cell_width = x_frac * (bbox[3] - bbox[1])
    y_frac = cell_side / (distance(float([bbox[1], bbox[2]]), float([bbox[1], bbox[4]]), units))
    cell_height = y_frac * (bbox[4] - bbox[2])

    xi = 0
    curr_x = bbox[1]

    while curr_x <= bbox[2]
        yi = 0
        curr_y = bbox[2]

        while curr_y <= bbox[4]
            cell1 = nothing
            cell2 = nothing

            if xi % 2 == 0 && yi % 2 == 0
                cell1 = Polygon([[[curr_x, curr_y], [curr_x, curr_y + cell_height],
                    [curr_x + cell_width, curr_y], [curr_x, curr_y]]])

                cell2 = Polygon([[[curr_x, curr_y + cell_height], [curr_x + cell_width, curr_y + cell_height],
                    [curr_x + cell_width, curr_y], [curr_x, curr_y + cell_height]]])
            elseif xi % 2 == 0 && yi % 2 == 1
                cell1 = Polygon([[[curr_x, curr_y], [curr_x + cell_width, curr_y + cell_height],
                    [curr_x + cell_width, curr_y], [curr_x, curr_y]]])

                cell2 = Polygon([[[curr_x, curr_y], [curr_x, curr_y + cell_height],
                    [curr_x + cell_width, curr_y + cell_height], [curr_x, curr_y]]])

            elseif yi % 2 == 0 && xi % 2 == 1
                cell1 = Polygon([[[curr_x, curr_y], [curr_x, curr_y + cell_height],
                    [curr_x + cell_width, curr_y + cell_height], [curr_x, curr_y]]])

                cell2 = Polygon([[[curr_x, curr_y], [curr_x + cell_width, curr_y + cell_height],
                    [curr_x + cell_width, curr_y], [curr_x, curr_y]]])
            elseif yi % 2 == 1 && xi % 2 == 1
                cell1 = Polygon([[[curr_x, curr_y], [curr_x, curr_y + cell_height],
                    [curr_x + cell_width, curr_y], [curr_x, curr_y]]])

                cell2 = Polygon([[[curr_x, curr_y + cell_height], [curr_x + cell_width, curr_y + cell_height],
                    [curr_x + cell_width, curr_y], [curr_x, curr_y + cell_height]]])
            end

            if mask != nothing
                # intersects...
                push!(results, cell1)
                push!(results, cell2)
            else
                push!(results, cell1)
                push!(results, cell2)
            end
            curr_y += cell_height
            yi += 1
        end
        xi += 1
        curr_x += cell_width
    end
    return FeatureCollection([Feature(x) for x in results])
end
