# Different Grids
This example showcases different types of grids on a map, delimited by a bbox (bounding box).

```
using Turf

box = [-95, 30 ,-85, 40]
cellside = 60

pg = point_grid(box, cellside) # a grid of points
rg = rectangle_grid(box, cellside, cellside / 2) # a grid of rectangles
tg = triangle_grid(box, cellside) # a grid of triangles
sg = square_grid(box, cellside) # a grid of squares

fc = [pg, rg, tg, sg]

for (i, fcoll) in enumerate(fc)
    open(joinpath(dirname(pathof(Turf)), "..", "docs") * "/examples/grids/grid$(i).result.geojson", "w") do file
        write(file, geojson(fcoll))
    end
end

```

![points](https://user-images.githubusercontent.com/40722053/61557883-bbc24280-aa65-11e9-9ae8-0db256a53828.JPG)     

![rect](https://user-images.githubusercontent.com/40722053/61557889-c4b31400-aa65-11e9-9541-1ecbbc7a024e.JPG)

![square](https://user-images.githubusercontent.com/40722053/61557902-d1376c80-aa65-11e9-9096-ee52df7254c4.JPG)     

![tri](https://user-images.githubusercontent.com/40722053/61557911-d85e7a80-aa65-11e9-9148-7e207906928a.JPG)
