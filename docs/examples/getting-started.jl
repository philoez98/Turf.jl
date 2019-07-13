using Turf

fc = GeoJSON.parsefile(pwd() * "/docs/examples/example.geojson")

centroid_point = centroid(fc)
center_point = center(fc)

push!(fc.features, Feature(centroid_point, Dict("marker-color" => "#fff000")))
push!(fc.features, Feature(center_point, Dict("marker-color" => "#eaa000")))

result = GeoJSON.geojson(fc)

open(pwd() * "/docs/examples/example.result.geojson", "w") do file
    write(file, result)
end
