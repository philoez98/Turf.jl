module Constants

export earthRadius, factors, unitsFactor, areaFactors

"""
Earth radius in meters (m).
"""
const earthRadius = 6371008.8


const factors = Dict(
    "centimeters" => earthRadius * 100,
    "centimetres" => earthRadius * 100,
    "degrees" => earthRadius / 111325,
    "feet" => earthRadius * 3.28084,
    "inches" => earthRadius * 39.370,
    "kilometers" => earthRadius / 1000,
    "kilometres" => earthRadius / 1000,
    "meters" => earthRadius,
    "metres" => earthRadius,
    "miles" => earthRadius / 1609.344,
    "millimeters" => earthRadius * 1000,
    "millimetres" => earthRadius * 1000,
    "nauticalmiles" => earthRadius / 1852,
    "radians" => 1,
    "yards" => earthRadius / 1.0936
)



"""
Units of measurement factors based on 1 meter.
"""
const unitsFactor = Dict(
    "centimeters" => 100,
    "centimetres" => 100,
    "degrees" => 1 / 111325,
    "feet" => 3.28084,
    "inches" => 39.370,
    "kilometers" => 1 / 1000,
    "kilometres" => 1 / 1000,
    "meters" => 1,
    "metres" => 1,
    "miles" => 1 / 1609.344,
    "millimeters" => 1000,
    "millimetres" => 1000,
    "nauticalmiles" => 1 / 1852,
    "radians" => 1 / earthRadius,
    "yards" => 1 / 1.0936
)

"""
Area of measurement factors based on 1 square meter.
"""
const areaFactors = Dict(
    "acres" => 0.000247105,
    "centimeters" => 10000,
    "centimetres" => 10000,
    "feet" => 10.763910417,
    "inches" => 1550.003100006,
    "kilometers" => 0.000001,
    "kilometres" => 0.000001,
    "meters" => 1,
    "metres" => 1,
    "miles" => 3.86e-7,
    "millimeters" => 1000000,
    "millimetres" => 1000000,
    "yards" => 1.195990046
)

end