
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
function do_map(x::Vector{T}, mask::Raster) where T <: Number
    missingval = T == Bool ? false : NaN
    map = fill(missingval, dims(mask); missingval)
    map[mask] .= x
    map
end

do_map(x, env::Environment) = do_map(x, env.mask)

get_climate(speciesmask::Raster, env::Environment) = (first(env.pca_maps)[speciesmask], last(env.pca_maps)[speciesmask])
get_climate(species::String, spec::Species, env::Environment) = get_climate(spec.ranges[At(species)], env)

function points_to_geo(xs, ys)
    length(xs) < 1 && return GI.MultiPoint([(0,0)])
    GI.MultiPoint(collect(zip(xs, ys)))
end

points_to_geo(points) = GI.MultiPoint(points)

# find the centroid in climate space of a species
function get_centroid(speciesmask, env::Environment)
    x, y = get_climate(speciesmask, env)
    mean(x), mean(y)
end

function find_range_shapes(spec::Species, env::Environment)
    cors = Float64[]
    xrange = Float64[]
    yrange = Float64[]
    for sp in spec.names
        x, y = get_climate(sp, spec, env)
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
gethull(species::String, spec::Species, env::Environment) = convexhull(points_to_geo(get_climate(species, spec, env)...))
hullarea(species::String,spec::Species, env::Environment) = LibGEOS.area(gethull(species, spec, env))

# Get the geographical centroid of a species
function geocentroids(name, spec::Species)
    a = spec.ranges[At(name)]
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

function makeweights(xs, ys, binsize = 0.1)
    makebins(i) = floor(minimum(i)-binsize):binsize:ceil(maximum(i)+binsize)
    points_in_cell(x, y, hist) = hist.weights[findfirst(>(x), hist.edges[1])-1, findfirst(>(y), hist.edges[2])-1]

    hist = fit(Histogram, (xs, ys), (makebins(xs), makebins(ys)))
    [1/points_in_cell(xs[i], ys[i], hist) for i in eachindex(xs, ys)]
end

function fitellipse(speciesname::String, spec::Species, env::Environment, sigma = 2; weightmap = nothing) 
    xs, ys = get_climate(speciesname, spec, env)
    length(xs) < 3 && return Ellipse(0, 0, 0, 0, 0)
    fit(Ellipse, xs, ys, sigma; weight = !isnothing(weightmap) ? weightmap[spec.ranges[At(speciesname)]] : ones(length(xs)))
end

# pick a random real grid cells as ellipse centers, rather than a random point
function pick_ellipse_center(env::Environment; on_real_point = false)    
    if on_real_point 
        centerpoint = rand(1:length(env.pca1))
        return(env.pca1[centerpoint], env.pca2[centerpoint])
    end
    rescale(rand(), first(env.bbox)...), rescale(rand(), last(env.bbox)...)
end

function sample_ellipse(harea = 1, env::Environment = env; max_iter = 1000, on_real_point = false)
    el = Ellipse(0,0,0,0,0)
    harea == 0 && return el
    ovrlp = 0
    failsafe = 0
    while ovrlp < 0.8 && (failsafe += 1) < max_iter
        el = rand(Ellipse, pick_ellipse_center(env; on_real_point)..., area = harea)
        ovrlp = overlap(el, env.chull) 
    end
    failsafe == max_iter && @show harea
    el
end