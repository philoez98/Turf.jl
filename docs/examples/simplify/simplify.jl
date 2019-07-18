using Turf

shape = GeoJSON.parsefile(joinpath(dirname(pathof(Turf)), "..", "docs") * "/examples/simplify/geometry.geojson")

simplify!(shape, 0.4)

result = GeoJSON.geojson(shape)

open(joinpath(dirname(pathof(Turf)), "..", "docs") * "/examples/simplify/geometry.result.geojson", "w") do file
    write(file, result)
end
