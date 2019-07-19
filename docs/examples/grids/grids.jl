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
