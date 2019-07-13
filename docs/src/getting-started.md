# Getting Started

This is a very short and simple example to get you started with Turf. It shows a little script that extracts some informations from a FeatureCollection and updates its content.

Suppose we have a bunch of points scattered around on a map. We want to know where the center and the centroid of all those points are.
In this case we'll use a `.geojson` file prepared with some random data. You can find all the files needed for this example inside the [examples](https://github.com/philoez98/Turf.jl/tree/master/docs/examples) folder.

So let's start by importing `Turf`. Note that `Turf` also exports the `GeoJSON.jl` and the `GeoInterface.jl` packages, so we don't need to import them:

```
using Turf
```

To read a `.geojson` file into a `String` we use the function `parsefile()` from the `GeoJSON.jl` package:

```
fc = GeoJSON.parsefile(pwd() * "/docs/examples/example.geojson")
```
So now we have a FeatureCollection object created from our file. We start to work on it by calculating both the center and the centroid of its features:

```
centroid_point = centroid(fc)
center_point = center(fc)
```

Both functions return a `Point`. Let's push those points in the original FeatureCollection so that we can later visualize them.
We also add the `color-marker` property to highlight them once plotted:

```
push!(fc.features, Feature(centroid_point, Dict("marker-color" => "#fff000"))) # yellow
push!(fc.features, Feature(center_point, Dict("marker-color" => "#eaa000"))) # orange
```
Now let's convert our FeatureCollection back to a `String` using the `geojson()` function:

```
result = GeoJSON.geojson(fc)
```

As a last step, we write the modified FeatureCollection into a new file, called `example.result.geojson`:

```
open(pwd() * "/docs/examples/example.result.geojson", "w") do file
    write(file, result)
end
```

And we're done! If we then plot the result here's how it looks like:

![result](https://user-images.githubusercontent.com/40722053/61170999-a171fb80-a571-11e9-8efe-bda1b98d2fa1.JPG)

You can open the `example.result.geojson` file in [geojson.io](http://geojson.io/) and see the same result as above.

Here's the entire script:

```
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
```


### Further Reading

For more examples and tutorials take a look at the [Examples](@ref) section.  
For a complete overview of all the methods available in Turf, take a look at the [Methods](@ref) section.
