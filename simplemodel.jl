
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

# A basic Ellipse struct
struct Ellipse
    center_x::Float64
    center_y::Float64
    length::Float64
    width::Float64
    angle::Float64 # Given in radians
end

# Let's load the environment data
datadir = "/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data"
ENV["RASTERDATASOURCES_PATH"] = joinpath(datadir, "Rasterdatasources")
# datadir = "/home/raf/Data/Biodiversity/Distributions"

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

# create maps to visualize the PCAs. These are kept in the analysis going forward

map1 = fill(NaN, dims(bioclim_sa); missingval=NaN)
map2 = fill(NaN, dims(bioclim_sa); missingval=NaN)
map1[sa_mask] .= pca1
map2[sa_mask] .= pca2
pcas = RasterStack((pca1=map1, pca2=map2))
Plots.plot(pcas)

# Load the bird shapefiles and pick the ones in South America

const GI = GeoInterface

shapefiles = [
    joinpath(datadir, "Birds/batch_1.shp"),
    joinpath(datadir, "Birds/batch_2.shp"),
    joinpath(datadir, "Birds/batch_3.shp"),
    joinpath(datadir, "Birds/batch_4.shp"),
    joinpath(datadir, "Birds/batch_5.shp"),
]

# sorry, this takes time
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


# Plot the diversity in climate space
f = Figure()
a, s = Makie.scatter(f[1,1], collect(zip(pca1, pca2)); markersize = 0.1, color = diversity[sa_mask], colormap = cgrad(:Spectral, rev = true))
Colorbar(f[1,2],s)
display(f)

# a function to rasterize a species by name
function get_speciesmask(name; geoms = sa_geoms, mask = sa_mask)
    ret = reduce(.|, map(findall(==(name), geoms.sci_name)) do i
        boolmask(geoms.geometry[i]; to = mask, boundary = :touches) .& mask
    end)
    rebuild(ret; name = name)
end

# Some testing plots

# names of all species
allspecies = unique(sa_geoms.sci_name)

# This plots the distribution of a random species on the map
Plots.plot(bioclim_sa.bio1, fill = :grey)
Plots.plot!(get_speciesmask(rand(allspecies)))

# Creates two side-by-side plots, one in geographic space, the other in climate space
# Shows the occurrence points of all the species defined by speciesnames
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
  #  leg = Makie.Legend(f[1,3], a, framevisible = false, colorrange = (1,10)) # TODO: there was a problem with the legend code
   # for i in 1:length(speciesnames)
    #    leg.entrygroups[][1][2][i].elements[1].markersize = 20
   # end
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

# This figure shows the geographic range size of ranges on the environmental centroids (in environment spaecs)
p = Figure()
a = Axis(p[1,1], aspect = DataAspect())
Makie.scatter!(a, collect(zip(pca1, pca2)); markersize = 0.1, color=(:grey, 0.5))
Makie.scatter!(a, allcentroids; markersize = 0.6, color = ranges, label = "rangesizes")
Makie.Colorbar(p[1,2])
p


# These functions create convex hulls in environmental space and calculate their area
function get_climatepoints(species)
    smask = get_speciesmask(species)
    xs, ys = get_climate(smask)
    length(xs) < 1 && return GI.MultiPoint([(0,0)])
    GI.MultiPoint(collect(zip(xs, ys)))
end

gethull(species) = convexhull(get_climatepoints(species))
hullarea(species) = LibGEOS.area(gethull(species))

# get hull area for all species
allhulls = hullarea.(allspecies)

# show the mean size of the climatic range in environmeltal space
# TODO: make this function general
p = Figure()
a = Axis(p[1,1], aspect = DataAspect())
Makie.scatter!(a, collect(zip(pca1, pca2)); markersize = 0.3, color=(:grey, 0.5))
Makie.scatter!(a, allcentroids; markersize = 4, color = allhulls, label = "rangesizes", colormap = Reverse(:Spectral))
Makie.Colorbar(p[1,2], colormap = Reverse(:Spectral))
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
    Makie.scatter!(axs[i], gc[inds], markersize = 2)
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

# els = bioclim_sa.bio1 .> 25
# Polygonize the masks
# This needs the improve_polygonize branch of GeometryOps
centers = map(lookup(els)) do lookup
    parent(Rasters.shiftlocus(Rasters.Center(), lookup))
end
polygons = GeometryOps.polygonize(centers..., els)
Rasters.rplot(els)
Makie.poly!(polygons)

for p in polygons
    Makie.poly!(p)
end

# Check the polygon areas - potential range sizes
res = map(abs ∘ step, centers)
min_viable_area = 20 # Maybe this is really small?
too_small_removed = filter(polygons) do polygon
    count(boolmask(polygon; res)) > min_pixels
end 
Rasters.rplot(els)
Makie.poly!(too_small_removed)

# Join areas withing some allowed dispersal radius
dilate(a) = Stencils.mapstencil(Stencils.Circle{4}(), a) do hood
    any(hood)
end
dilated_polygons = GeometryOps.polygonize(centers..., dilate(els))
Rasters.rplot(els)
Makie.poly!(dilated_polygons)

multi_polygons = GI.MultiPolygon.(map(dilated_polygons) do dp
    filter(polygons) do p
        LibGEOS.within(p, dp)
    end 
end |> filter(x -> length(x) > 0));

# Plot habitat groups that may be the same species
Makie.poly(dilated_polygons)
for mp in multi_polygons
    Makie.poly!(mp)
end

p = Rasters.rplot(pcas.pca1; transparency=true)
Makie.poly!.(p.axis, multi_polygons)
Makie.save("species.png", p)
p = Makie.scatter( collect(zip(pca1, pca2)); markersize = 0.1, color = :grey)
collect(zip(pcas.pca1, pcas.pca2))[els]
Makie.scatter!(collect(zip(pcas.pca1, pcas.pca2))[els]; markersize = 0.1, color = :red)
Makie.save("pca.png", p)


# Now rasterise these into the PCAs? 
