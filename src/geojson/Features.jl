module Features

include("Geometries.jl")

using .Geometries

export Feature, FeatureCollection

abstract type AbstractFeature end

const Properties = Union{Dict{String, Any}, Nothing}

mutable struct Feature <: AbstractFeature
        type::String
        geometry::Geometry
        properties::Properties

        Feature(type::String, geometry::Geometry, properties::Properties) = type === "Feature" ? new(type, geometry, properties) : setType(new(type, geometry, properties))
end

Feature(type::String, geometry::Geometry) = Feature(type, geometry, nothing)
setType(feat::Feature) = feat.type = "Feature"

mutable struct FeatureCollection <: AbstractFeature
        type::String
        features::Union{Vector{Feature}}

        FeatureCollection(type::String, features::Vector{Feature}) = type === "FeatureCollection" ? new(type, features) : setType(new(type, features))
end
FeatureCollection() = FeatureCollection("FeatureCollection", [])
setType(feat::FeatureCollection) = feat.type = "FeatureCollection"

end # module
