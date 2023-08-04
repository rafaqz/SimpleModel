using Rasters, RasterDataSources
using Stencils
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

struct Ellipse
    center_x::Float64
    center_y::Float64
    length::Float64
    width::Float64
    angle::Float64 # Given in radians
end

# datadir = "/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data"
# ENV["RASTERDATASOURCES_PATH"] = joinpath(datadir, "Rasterdatasources")
datadir = "/home/raf/Data/Biodiversity/Distributions"

# download the bioclim variables for south america

# bioclim = RasterStack(CHELSA{BioClim}; lazy=true)
# bioclim_sa = bioclim[X=-89 .. -33, Y=-57 .. 13]
# bioclim_sa = Rasters.aggregate(mean, replace_missing(bioclim_sa, NaN), 10)

# Use WorldClim for now for speed
bioclim = RasterStack(WorldClim{BioClim}; lazy=true)
bioclim_sa = bioclim[X=-89 .. -33, Y=-57 .. 13]
# bioclim_sa = Rasters.aggregate(mean, replace_missing(bioclim_sa, NaN), 10)
sa_mask = boolmask(bioclim_sa.bio15)
Plots.plot(bioclim_sa.bio1)

# Fit a PCA model to the climate and extract the two primary components

# logged = RasterStack(Tuple(map(x -> x[1] > 11 ? log.(x[2]) : x[2], enumerate(values(replace_missing(bioclim_sa))))))
big_mat = reduce(vcat, map(enumerate(bioclim_sa)) do (i, A)
    permutedims(zscore(i > 11 ? log.(A[sa_mask] .+ 1) : A[sa_mask]))
end)
# histogram(big_mat'; bins=20, ticks=false, label=false, title=(1:19)')

model = fit(PCA, big_mat; maxoutdim=2)
pred = MultivariateStats.transform(model, big_mat)
pca1 = pred[1, :]
pca2 = pred[2, :]

# visualize the PCA biplot

p = Makie.scatter(collect(zip(pca1, pca2)); color=(:grey, 0.5))
for i in 1:19
    Makie.lines!([(0,0), 15 .* (projection(model)[i,:]...,)], color = :red)
    Makie.text!((15 .* projection(model)[i,:])..., text = string(names(bioclim)[i]))
end
p

# create maps to visualize the PCAs

map1 = fill(NaN, dims(bioclim_sa); missingval=NaN)
map2 = fill(NaN, dims(bioclim_sa); missingval=NaN)
map1[sa_mask] .= pca1
map2[sa_mask] .= pca2
pcas = RasterStack((pca1=map1, pca2=map2))
Plots.plot(pcas)

# Load the bird shapefiles and pick the ones in South America

using Shapefile, GeoInterface, Extents#, ProgressMeter
const GI = GeoInterface
using GeometryOps

shapefiles = [
    joinpath(datadir, "Birds/batch_1.shp"),
    joinpath(datadir, "Birds/batch_2.shp"),
    joinpath(datadir, "Birds/batch_3.shp"),
    joinpath(datadir, "Birds/batch_4.shp"),
    joinpath(datadir, "Birds/batch_5.shp"),
]

sa_geoms = reduce(vcat, map(shapefiles) do sf
    df = DataFrame(Shapefile.Table(sf))
    filter(df) do row
        ext = GI.calc_extent(GI.trait(row.geometry), row.geometry)
        Extents.intersects(ext, Extents.extent(bioclim_sa))
    end
end)

# Count total diversity
diversity = rasterize(count, sa_geoms; to=sa_mask, boundary=:touches)
diversity .*= sa_mask
Plots.plot(diversity)
Plots.plot(diversity; clims=(400, 650))


# Plot the diversity in climate space
f = Figure()
a, s = Makie.scatter(f[1,1], collect(zip(pca1, pca2)); markersize = 0.1, color = diversity[sa_mask], colormap = cgrad(:Spectral, rev = true))
Colorbar(f[1,2],s)
display(f)

# a function to rasterize a species by name

function get_speciesmask(name; geoms = sa_geoms, mask = sa_mask)
    reduce(.|, map(findall(==(name), geoms.sci_name)) do i
        boolmask(geoms.geometry[i]; to = mask, boundary = :touches) .& mask
    end)
end

# Some testing plots
allspecies = unique(sa_geoms.sci_name)
Plots.plot(get_speciesmask(rand(allspecies)))

function plot_species(speciesnames)
    f = Figure(resolution = (1500, 700))
    a = Axis(f[1,1], aspect = DataAspect())
    b = Axis(f[1,2], aspect = DataAspect())
    Makie.plot!(a, sa_mask, colormap = "Greys")
    Makie.scatter!(b, collect(zip(pca1, pca2)); markersize = 0.1, color=(:grey, 0.5))
    colsize!(f.layout, 1, Auto(0.73))
    rich = fill(NaN, dims(sa_mask))

    for (i, spec) in enumerate(speciesnames)
        species_mask = get_speciesmask(spec)
        rich[species_mask] .= i
        Makie.scatter!(b, map1[species_mask], map2[species_mask]; 
            markersize = 0.2, label = spec,
            colormap=(:tab10, 0.5), color=i, colorrange=(1, 10)
        )
    end

    Makie.plot!(a, rich, colormap=(:tab10, 0.5), colorrange=(1, 10))
    leg = Makie.Legend(f[1,3], b, framevisible = false)
    for i in 1:length(speciesnames)
        leg.entrygroups[][1][2][i].elements[1].markersize = 20
    end
    f
end

# plot some random species and look at their distribution
specs = rand(allspecies, 5)
plot_species(specs)

function get_climate(speciesmask; pcas = (map1, map2))
    filter(isfinite, map1[speciesmask]), filter(isfinite, map2[speciesmask])
end

# find the centroid in climate space of a species
function get_centroid(speciesmask; pcas = (map1, map2))
    x, y = get_climate(speciesmask; pcas)
    mean(x), mean(y)
end

# plot the density of centroids
allcentroids = get_centroid.(get_speciesmask.(allspecies))
ranges = count.(get_speciesmask.(allspecies))

p = Figure()
a = Axis(p[1,1], aspect = DataAspect())
Makie.scatter!(a, collect(zip(pca1, pca2)); markersize = 0.1, color=(:grey, 0.5))
Makie.scatter!(a, allcentroids; markersize = 0.2, color = ranges, label = "rangesizes")
Makie.Colorbar(p[1,2])
p



# Get the geographical centroid of a species
function geocentroids(name)
    a = get_speciesmask(name)
    dp = DimPoints(a)[a]
    isempty(dp) && return(NaN, NaN)
    mean(first, dp), mean(last, dp)
end

gc = geocentroids.(allspecies)


# plot the location of all the centroids on the map
f = Figure()
a = Axis(f[1,1], aspect = DataAspect())
Makie.plot!(a, sa_mask, colormap = :Greys)
Makie.scatter!(gc, color = ranges)
f

# divide the rangesizes into quartiles for plotting
function asquantile(vec, n)
    q = quantile(vec, (0:n-1)./n)
    [findlast(<(r), q) for r in vec]
end

rangequants = asquantile(ranges, 4)

# and plot them

f = Figure()
aspect = DataAspect()
axs = [Axis(f[1,1]; aspect), Axis(f[2,1]; aspect), Axis(f[1,2]; aspect), Axis(f[2,2]; aspect)]
for i in 1:4
    inds = findall(==(i), rangequants)
    Makie.plot!(axs[i], sa_mask, colormap = :Greys)
    Makie.scatter!(axs[i], gc[inds], markersize = 0.9)
end
f

# save("climate_space_5_birds.png")

background = fit(Histogram, (pca1, pca2),  (-10:0.3:10, -10:0.3:10))
spec = rand(allspecies)
x,y = get_climate(get_speciesmask(spec))
sp = fit(Histogram, (x,y),  (-10:0.3:10, -10:0.3:10))
Plots.plot(Plots.plot(sp, aspect_ratio = 1), Plots.plot(background, aspect_ratio = 1), layout = 2)
rel = sp.weights ./ background.weights
replace!(rel, 0 => NaN)
Plots.heatmap(rel')




# not really necessary
#using GeometryBasics
#simplified = GeometryOps.simplify(sa_geoms[1:100]; tol=0.0000001);
#GeoInterface.convert(GeometryBasics, simplified)
#rast = rasterize(count, sa_geoms; to=sa_mask, threaded=false)
#Plots.plot(rast)

# check if the point `(xp, yp)`` is in an ellipse centered on (`x`,`y`) 
# with angle `an` (in radians), width `a` and height `b`
function in_ellipse((xp,yp), an, x, y, a, b)
    a = (cos(an) * (xp-x) + sin(an)*(yp-y))^2 / a^2
    b = (sin(an) * (xp-x) + cos(an)*(yp-y))^2 / b^2
    a+b < 1
end

function distance(point::Tuple, el::Ellipse)
    cosa = cos(el.angle)
    sina = sin(el.angle)
    rel_x = first(point) - el.center_x
    rel_y = last(point) - el.center_y
    a = (cosa * rel_x + sina * rel_y)^2 / el.length^2
    b = (sina * rel_x + cosa * rel_y)^2 / el.width^2
    a+b
end

in_ellipse(point, el::Ellipse) = distance(point, el) <= 1

randlims(lims) = (last(lims)-first(lims))*rand()+first(lims)

import Random.rand
function Random.rand(::Type{Ellipse}; 
    xlims=(0,1), ylims=(0,1), lengthlims=(0.01,1), widthlims=(0.01,1)
) 
    Ellipse(randlims(xlims), randlims(ylims), randlims(lengthlims), randlims(widthlims), rand()π)
end

# Test the ellipse code
using Plots
pts = vec(collect(Iterators.product(1:5:1000, 1:5:1000)))
el = rand(Ellipse; xlims=(1, 1000), ylims=(1, 1000), lengthlims=(1,400), widthlims=(1,400))
Plots.scatter(pts; aspect_ratio = 1, marker_z = [in_ellipse(x, el) for x in pts], msw = 0, ms = 1)

el = rand(Ellipse; xlims = (-5,10), ylims = (-5,5), lengthlims = 1.5, widthlims = 1.5)
els = [in_ellipse(pt, el) for pt in zip(pcas.pca1, pcas.pca2)]
h = Makie.scatter(collect(zip(pca1, pca2)); markersize = 0.1, color = :grey)
collect(zip(pcas.pca1, pcas.pca2))[els]
Makie.scatter!(collect(zip(pcas.pca1, pcas.pca2))[els]; markersize = 0.1, color = :red)

# Project niche pca elipses back to real space as masks
Rasters.rplot(els)

# Polygonize the masks
# This needs the improve_polygonize branch of GeometryOps
using GeometryOps, GeometryBasics, GLMakie, GeoInterface
# This is still buggy...
centers = map(lookup(els)) do lookup
    parent(Rasters.shiftlocus(Rasters.Center(), lookup))
end
@time polygons = GeometryOps.polygonize(centers..., els);
Rasters.rplot(els)
# Slightly buggy still, some holes are missed 
Makie.poly!(polygons)

# Check the polygon areas - potential range sizes
areas = map(polygons) do polygon
    res = map(abs ∘ step, centers)
    count(boolmask(polygon; res)) * prod(res) 
end |> sort
# Sort the areas 

# Now rasterise these into the PCAs? 
