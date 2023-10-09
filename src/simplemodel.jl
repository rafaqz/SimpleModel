
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

# Create a raster map from a vector of observations in the right order
function do_map(pca, sa_mask; missingval = NaN)
    map = fill(missingval, dims(sa_mask); missingval)
    map[sa_mask] .= pca
    map
end



# a function to rasterize a species by name
function get_speciesmask(name; geoms = sa_geoms, mask = sa_mask)
    ret = reduce(.|, map(findall(==(name), geoms.sci_name)) do i
        boolmask(geoms.geometry[i]; to = mask, boundary = :touches) .& mask
    end)
    rebuild(ret; name = name)
end

function get_climate(speciesmask::Raster; pcas = pca_maps)
    filter(isfinite, first(pcas)[speciesmask]), filter(isfinite, last(pcas)[speciesmask])
end

get_climate(species::String; pcas = pca_maps, allranges = allranges) = get_climate(allranges[At(species)]; pcas)

function points_to_geo(xs, ys)
    length(xs) < 1 && return GI.MultiPoint([(0,0)])
    GI.MultiPoint(collect(zip(xs, ys)))
end

points_to_geo(points) = GI.MultiPoint(points)

# find the centroid in climate space of a species
function get_centroid(speciesmask; pcas = pca_maps)
    x, y = get_climate(speciesmask; pcas)
    mean(x), mean(y)
end

# These functions create convex hulls in environmental space and calculate their area
gethull(species) = convexhull(points_to_geo(get_climate(species)...))
hullarea(species) = LibGEOS.area(gethull(species))

# Get the geographical centroid of a species
function geocentroids(name)
    a = allranges[At(name)]
    dp = DimPoints(a)[a]
    isempty(dp) && return(NaN, NaN)
    mean(first, dp), mean(last, dp)
end

# divide a vector into quartiles
function asquantile(vec, n)
    q = quantile(vec, (0:n-1)./n)
    [findlast(<(r), q) for r in vec]
end

function overlap(el::Ellipse, polygon; n = 100)
    ellipse_points = points_to_geo(coordinates(el, n))
    ellipse_poly = GI.Polygon([GI.LinearRing(GI.getpoint(ellipse_points))])
    LibGEOS.area(LibGEOS.intersection(polygon, ellipse_poly)) / LibGEOS.area(ellipse_poly)
end

points_in_cell(x, y) = hist2.weights[findfirst(>(x), hist2.edges[1])-1, findfirst(>(y), hist2.edges[2])-1]
getweights(xs, ys) = [1/points_in_cell(xs[i], ys[i]) for i in eachindex(xs, ys)]
function fitellipse(speciesname::String, sigma = 2; weighted = false) 
    xs, ys = get_climate(speciesname)
    length(xs) < 3 && return Ellipse(0, 0, 0, 0, 0)
    fit(Ellipse, xs, ys, sigma; weight = weighted ? getweights(xs, ys) : ones(length(xs)))
end