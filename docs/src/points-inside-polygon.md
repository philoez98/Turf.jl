# Points inside a Polygon

Example that shows how to find all points that are inside a specific shape (polygon) on the map.


```

using Turf

geo_data = GeoJSON.parsefile(joinpath(dirname(pathof(Turf)), "..", "docs") * "/examples/points-inside-polygon/points-inside-polygon.geojson")

polygon = geo_data.features[1].geometry

for i in 2:length(geo_data.features)
    feature = geo_data.features[i]
    point = feature.geometry

    if within(point, polygon)
        feature.properties = Dict("marker-color" => "#ff0000")
    end
end

result = geojson(geo_data)

open(joinpath(dirname(pathof(Turf)), "..", "docs") * "/examples/points-inside-polygon/points-inside-polygon.result.geojson", "w") do file
    write(file, result)
end

```

![inside](https://user-images.githubusercontent.com/40722053/61465345-2f862180-a978-11e9-9b52-9a2488ebc1c3.JPG)
