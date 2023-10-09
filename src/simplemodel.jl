
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
function do_map(x, mask; missingval = NaN)
    map = fill(missingval, dims(mask); missingval)
    map[mask] .= x
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

function find_range_shapes(allspecies)
    cors = Float64[]
    xrange = Float64[]
    yrange = Float64[]
    for spec in allspecies
        x, y = get_climate(spec)
        if isempty(x)
            push!(cors, 0)
            push!(xrange, 0)
            push!(yrange, 0)
        else
            push!(cors, cor(x,y))
            push!(xrange, maximum(x) - minimum(x))
            push!(yrange, maximum(y) - minimum(y))
        end
    end
    cors, xrange, yrange
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

function makeweights(x, y, binsize = 0.1)
    makebins(i) = floor(minimum(i)):binsize:ceil(maximum(i))
    points_in_cell(x, y, hist) = hist.weights[findfirst(>(x), hist.edges[1])-1, findfirst(>(y), hist.edges[2])-1]

    hist = fit(Histogram, (x, y), (makebins(x), makebins(y)))
    [1/points_in_cell(x[i], y[i], hist) for i in eachindex(x, y)]
end

function fitellipse(speciesname::String, sigma = 2; weighted = false) 
    xs, ys = get_climate(speciesname)
    length(xs) < 3 && return Ellipse(0, 0, 0, 0, 0)
    fit(Ellipse, xs, ys, sigma; weight = weighted ? weightmap[allranges[speciesname]] : ones(length(xs)))
end

# pick a random real grid cells as ellipse centers, rather than a random point
function pick_ellipse_center(x = pca1, y = pca2, bbox_all = bbox_all; on_real_point = false)    
    if on_real_point 
        centerpoint = rand(1:length(x))
        return(x[centerpoint], y[centerpoint])
    end
    rescale(rand(), first(bbox_all)...), rescale(rand(), last(bbox_all)...)
end

function sample_ellipse(harea = 1; max_iter = 1000, on_real_point = false)
    el = Ellipse(0,0,0,0,0)
    harea == 0 && return el
    ovrlp = 0
    failsafe = 0
    while ovrlp < 0.8 && (failsafe += 1) < max_iter
        el = rand(Ellipse, pick_ellipse_center(;on_real_point)..., area = harea)
        ovrlp = overlap(el, chull_all) 
    end
    failsafe == max_iter && @show harea
    el
end