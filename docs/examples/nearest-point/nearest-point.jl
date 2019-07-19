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
