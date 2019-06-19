module GeoJSON

include("Features.jl")

using .Features

export GeoJson

mutable struct GeoJson
    content::Union{FeatureCollection, Feature}
end

GeoJson() = GeoJson(FeatureCollection())


end # module
