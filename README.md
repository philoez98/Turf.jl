# Turf.jl

| **Documentation** | **Build Status** |
|:---:|:---:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://philoez98.github.io/Turf.jl/stable) | [![Linux/MacOS](https://travis-ci.org/philoez98/Turf.jl.svg?branch=master)](https://travis-ci.org/philoez98/Turf.jl)  [![Windows](https://ci.appveyor.com/api/projects/status/deghewsv2gra487s?svg=true)](https://ci.appveyor.com/project/philoez98/turf-jl)  [![Coverage Status](https://coveralls.io/repos/github/philoez98/Turf.jl/badge.svg?branch=master)](https://coveralls.io/github/philoez98/Turf.jl?branch=master)  [![codecov](https://codecov.io/gh/philoez98/Turf.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/philoez98/Turf.jl) |


A spatial analysis library written in Julia, ported from the great [*Turf.js*](https://github.com/Turfjs/turf).

Turf.jl uses [GeoInterface.jl](https://github.com/JuliaGeo/GeoInterface.jl) and [GeoJSON.jl](https://github.com/JuliaGeo/GeoJSON.jl) to create and handle all geographic data.


## Installation

Turf.jl can be installed either from the REPL using Pkg mode:

```
pkg> add Turf
```
or via `Pkg`:

```
julia> import Pkg; Pkg.add("Turf")
```
Alternatively if you want to use the latest version available, you can do:

```
pkg> add Turf#master
```

## Example

As an example, let's try to identify which points are within a certain polygon on the map and mark them with a different color.
We can do this using Turf.

```julia
# Turf already exports all symbols of GeoInterface.jl and GeoJSON.jl, so there's no need to import them
using Turf

# Let's create a basic geojson with a polygon and some points
fc = """
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "properties": {},
      "geometry": {
        "type": "Polygon",
        "coordinates": [
          [
            [5.185547, 47.753949],
            [-0.703125, 39.637989],
            [7.910156, 30.974436],
            [12.568359, 48.107339],
            [5.185547, 47.753949]
          ]
        ]
      }
    },
     {
      "type": "Feature",
      "properties": {},
      "geometry": {"type": "Point", "coordinates": [12.304688, 43.3882]}
    },
    {
      "type": "Feature",
      "properties": {},
      "geometry": {"type": "Point", "coordinates": [8.964844, 42.031854]}
    },
    {
      "type": "Feature",
      "properties": {},
      "geometry": {"type": "Point", "coordinates": [2.548828, 41.901134]}
    },
    {
      "type": "Feature",
      "properties": {},
      "geometry": {"type": "Point", "coordinates": [9.84375, 36.525174]}
    },
    {
      "type": "Feature",
      "properties": {},
      "geometry": {"type": "Point", "coordinates": [-0.878906, 38.546418]}
    },
    {
      "type": "Feature",
      "properties": {},
      "geometry": {"type": "Point", "coordinates": [3.47168, 46.407564]}
    }
  ]
}
"""
# convert the geojson into a FeatureCollection object
geo_data = GeoJSON.parse(fc)

# extract the polygon
polygon = geo_data.features[1].geometry

# loop to see which points are inside the polygon
# and add a 'marker-color' property to them
for i in 2:length(geo_data.features)
    feature = geo_data.features[i]
    point = feature.geometry

    if within(point, polygon)
        feature.properties = Dict("marker-color" => "#ff0000")
    end
end

# convert the FeatureCollection back to geojson
result = geojson(geo_data)

```

If we then plot the results here's what we get:

before:

![before (2)](https://user-images.githubusercontent.com/40722053/60754992-a4a53e80-9fe9-11e9-98d5-9bd889fcb0f0.JPG)

after:

![after](https://user-images.githubusercontent.com/40722053/60755010-e33af900-9fe9-11e9-89d9-2e3164e4a7ca.JPG)

For more examples, please take a look at the [Examples](https://philoez98.github.io/Turf.jl/latest/examples/) section of the documentation.


## Documentation

- [**STABLE**](https://philoez98.github.io/Turf.jl/stable): most recent tagged version.
- [**DEVEL**](https://philoez98.github.io/Turf.jl/latest): development version.


## Getting Help

- **Have a bug to report?** [Open an issue](https://github.com/philoez98/Turf.jl/issues/new/choose). Include the version of Turf and Julia, a full log, and some code that shows the issue.
- **Have a feature request?** [Open an issue](https://github.com/philoez98/Turf.jl/issues/new/choose). Tell us what the feature should do and why you want the feature.

## Available Functionality

A list with the currently available features can be found [here](https://github.com/philoez98/Turf.jl/blob/master/Turf.md).
Please open an issue if there's a specific *Turf.js* method that you'd like to have implemented.

## Contribute

Contributions are highly welcomed and appreciated.
If you want to contribute to this project, feel free to open an [issue](https://github.com/philoez98/Turf.jl/issues/new/choose) to discuss your proposal and its implementation. Once it's all good make a [pull request](https://github.com/philoez98/Turf.jl/pulls) and you're done.
Thank you for helping out!
