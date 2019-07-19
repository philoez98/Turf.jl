# Simplify a Geometry

This example shows how to simplify any geometry. Adjust the `tolerance` parameter to tune how much simpler the shape should be.

```

using Turf

shape = GeoJSON.parsefile(joinpath(dirname(pathof(Turf)), "..", "docs") * "/examples/simplify/geometry.geojson")

simplify!(shape, 0.4)

result = GeoJSON.geojson(shape)

open(joinpath(dirname(pathof(Turf)), "..", "docs") * "/examples/simplify/geometry.result.geojson", "w") do file
    write(file, result)
end

```

original:

![not_sim](https://user-images.githubusercontent.com/40722053/61466434-331aa800-a97a-11e9-8a00-da7b5b9d872c.JPG)

simplified:

![sim](https://user-images.githubusercontent.com/40722053/61466480-4463b480-a97a-11e9-9dae-498f2140467a.JPG)
