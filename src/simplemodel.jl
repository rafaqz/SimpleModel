
# Load a shitload of packages
using Rasters
using RasterDataSources
using MultivariateStats
using Distributions
using Plots
using Random
using StatsBase
using GLMakie
using Statistics
using DataFrames
using ArchGDAL
using GeoInterfaceMakie
using GeoInterfaceRecipes
using GeometryOps
using GeometryBasics
using GeoInterface
using LibGEOS
using Stencils
using Extents
using Shapefile

# a function to rasterize a species by name
function get_speciesmask(name; geoms = sa_geoms, mask = sa_mask)
    ret = reduce(.|, map(findall(==(name), geoms.sci_name)) do i
        boolmask(geoms.geometry[i]; to = mask, boundary = :touches) .& mask
    end)
    rebuild(ret; name = name)
end

function get_climate(speciesmask; pcas = pca_maps)
    filter(isfinite, first(pcas)[speciesmask]), filter(isfinite, last(pcas)[speciesmask])
end

get_climate(species; pcas = pca_maps) = get_climate(get_speciesmask(species); pcas)

function points_to_geo(xs, ys)
    length(xs) < 1 && return GI.MultiPoint([(0,0)])
    GI.MultiPoint(collect(zip(xs, ys)))
end

# find the centroid in climate space of a species
function get_centroid(speciesmask; pcas = (map1, map2))
    x, y = get_climate(speciesmask; pcas)
    mean(x), mean(y)
end

# These functions create convex hulls in environmental space and calculate their area
gethull(species) = convexhull(points_to_geo(get_climatepoints(species)...))
hullarea(species) = LibGEOS.area(gethull(species))

# Get the geographical centroid of a species
function geocentroids(name)
    a = get_speciesmask(name)
    dp = DimPoints(a)[a]
    isempty(dp) && return(NaN, NaN)
    mean(first, dp), mean(last, dp)
end

# divide a vector into quartiles
function asquantile(vec, n)
    q = quantile(vec, (0:n-1)./n)
    [findlast(<(r), q) for r in vec]
end
