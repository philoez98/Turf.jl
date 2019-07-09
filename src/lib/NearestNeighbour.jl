
struct NNStatistics
    units::String
    arealUnits::String
    area::Polygon
    observedMeanDistance::Real
    expectedMeanDistance::Real
    nearestNeighbourIndex::Real
    pointsCount::Real
    zScore::Real
end


"""
    analysis(data::F, area::Union{P, Nothing}=nothing, units::String="kilometers") where {F <: AbstractFeatureCollection, P <: AbstractPolygon}

Nearest Neighbor Analysis calculates an index based the average distances
between points in the dataset, thereby providing inference as to whether the
data is clustered, dispersed, or randomly distributed within the study area.

It returns a Polygon of the study area, with the results of
the analysis attached as part of of the `nearestNeighborAnalysis` property
of the study area's `properties`. The attached
[_z_-score](https://en.wikipedia.org/wiki/Standard_score) indicates how many
standard deviations above or below the expected mean distance the data's
observed mean distance is. The more negative, the more clustered. The more
positive, the more evenly dispersed. A _z_-score between -2 and 2 indicates
a seemingly random distribution. That is, within _p_ of less than 0.05, the
distribution appears statistically significantly neither clustered nor
dispersed.

**Remarks**

- Though the analysis will work on any FeatureCollection type, it
works best with Point collections.

- This analysis is **very** sensitive to the study area provided.
If no Polygon is passed as the study area, the function draws a box
around the data, which may distort the findings. This analysis works best
with a bounded area of interest within with the data is either clustered,
dispersed, or randomly distributed. For example, a city's subway stops may
look extremely clustered if the study area is an entire state. On the other
hand, they may look rather evenly dispersed if the study area is limited to
the city's downtown.
"""
function analysis(data::F, area::Union{P, Nothing}=nothing, units::String="kilometers") where {F <: AbstractFeatureCollection, P <: AbstractPolygon}
    area == nothing && (area = bbox_polygon(bbox(data)[1]))

    feats = []

    for ft in data.features
        push!(feats, centroid(ft.geometry))
    end

    n = length(feats)

    obsDist = []
    for (index, p) in enumerate(feats)
        other = filter((f, i) -> i !== index, feats)
        push!(obsDist, distance(p, nearestpoint(p, other).coordinates, units))
    end

    reduce((a, b) -> a + b, obsDist)

    popDensity = n / convert_area(area(area), "meters", units)
    expDist = 1 / (2 * sqrt(popDensity))
    var = 0.26136 / (sqrt(n * popDensity))

    return NNStatistics(units, units * "²", area, obsDist, expDist, obsDist / expDist, n, (obsDist - expDist) / var)

end
