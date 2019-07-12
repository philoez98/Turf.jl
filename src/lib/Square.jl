"""
    square(bbox::Vector{T}) where {T <: Real}

Take a bounding box and calculates the minimum square bounding box that would contain the input.

# Examples
```julia
julia> using Turf

julia> bbox = [-1, 1, 2, 3.5]
4-element Array{Float64,1}:
 -1.0
  1.0
  2.0
  3.5

julia> square(bbox)
4-element Array{Float64,1}:
 -1.0
  0.75
  2.0
  3.75
```
"""
function square(bbox::Vector{T}) where {T <: Real}
    west = bbox[1]
    south = bbox[2]
    east = bbox[3]
    north = bbox[4]

    hDist = distance(float(bbox[1:2]), float([east, south]))
    vDist = distance(float(bbox[1:2]), float([west, north]))

    if hDist >= vDist
        vMidPoint = (south + north) / 2

        return [west, vMidPoint - ((east - west) / 2), east, vMidPoint + ((east - west) / 2)]
    else
        hMidPoint = (west + east) / 2

        return [hMidPoint - ((north - south) / 2), south, hMidPoint + ((north - south) / 2),  north]
    end
end
