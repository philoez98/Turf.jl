# Nearest Point

This example shows how to find the nearest point to a specific location on the map. The destination is shown in red.

```
using Turf

fc = GeoJSON.parsefile(joinpath(dirname(pathof(Turf)), "..", "docs") * "/examples/nearest-point/points.geojson")

function get_point_with_property(features, prop)
    point = nothing
    for (i, feature) in enumerate(features)
        if haskey(feature.properties, prop)
            point = features[i]
            deleteat!(features, i)
        end
    end
    return point
end

destination = get_point_with_property(fc.features, "tag")

nearest = nearestpoint(destination.geometry, [x.geometry for x in fc.features])

index = findfirst(x -> x.geometry.coordinates == nearest.coordinates, fc.features)

fc.features[index].properties = Dict("marker-color" => "#fbb000", "tag" => "nearest-point") # yellow marker
push!(fc.features, destination)

result = geojson(fc)

open(joinpath(dirname(pathof(Turf)), "..", "docs") * "/examples/nearest-point/points.result.geojson", "w") do file
    write(file, result)
end

```

![near](https://user-images.githubusercontent.com/40722053/61554099-75b4b100-aa5c-11e9-909f-57f9d8f82969.JPG)
