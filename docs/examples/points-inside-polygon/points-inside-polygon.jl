using Turf

geo_data = GeoJSON.parsefile(pwd() * "/docs/examples/points-inside-polygon/points-inside-polygon.geojson")

polygon = geo_data.features[1].geometry

for i in 2:length(geo_data.features)
    feature = geo_data.features[i]
    point = feature.geometry

    if within(point, polygon)
        feature.properties = Dict("marker-color" => "#ff0000")
    end
end

result = geojson(geo_data)

open(pwd() * "/docs/examples/points-inside-polygon/points-inside-polygon.result.geojson", "w") do file
    write(file, result)
end
