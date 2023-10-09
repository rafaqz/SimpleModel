using Plots
using GLMakie
using Rasters
using JLSO
using GeoInterface; const GI = GeoInterface
using LibGEOS
using ConcaveHull


include("plotting.jl")
include("ellipse.jl")
include("simplemodel.jl")

###--- First we load all the data

datadir = "/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data"

obj = try
    JLSO.load(joinpath(datadir, "processed_objects.jls"))
catch
    include("prepare_data.jl")

    ## Get the environmental data
    bioclim_sa = prepare_environment(datadir)
    sa_mask = boolmask(bioclim_sa.bio15)

    # and visualize
    Plots.plot(bioclim_sa.bio1)

    ## get the PCA
    pca1, pca2, loads = do_pca(bioclim_sa, sa_mask)

    # and convert the results to raster
    pca_maps = RasterStack((pca1=do_map(pca1, sa_mask), pca2=do_map(pca2, sa_mask)))

    # and visualize
    Plots.plot(pca_maps)
    biplot(pca1, pca2, loads, string.(names(bioclim_sa)))

    # Get the bird data (this takes time)
    sa_geoms = loadranges("Birds", 5, sa_mask, "/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data")
    # names of all species
    allspecies = unique(sa_geoms.sci_name)
    allranges = RasterSeries([get_speciesmask(name; geoms = sa_geoms, mask = sa_mask) for name in allspecies], (; name = allspecies))
    obj = :obj => (pca1, pca2, pca_maps, bioclim_sa, sa_mask, allranges, allspecies)
    JLSO.save(joinpath(datadir, "processed_objects.jls"), obj)
    Dict(obj)
end

pca1, pca2, pca_maps, bioclim_sa, sa_mask, allranges, allspecies = obj[:obj]
###--- Exploratory data analysis

# Count total diversity
diversity = reduce(+, allranges)
Plots.plot(diversity)


# Plot the diversity in climate space
f = Figure()
a, s = Makie.scatter(f[1,1], collect(zip(pca1, pca2)); markersize = 0.1, color = diversity[sa_mask], colormap = cgrad(:Spectral, rev = true))
Colorbar(f[1,2],s)
display(f)


# plot some random species and look at their distribution
specs = rand(allspecies, 5)
plot_species(specs)


## how are ranges shaped?
cors, xrange, yrange = find_range_shapes(allspecies)
histogram(cors)
Plots.scatter(xrange, yrange)


###--- Patterns of range and niche size in climate space

# get environmental centroids and range sizes for all species
allcentroids = vec(get_centroid.(allranges))
ranges = vec(count.(allranges))

# This figure shows the geographic range size of ranges on the environmental centroids (in environment spaecs)
overplot_pca_space(allcentroids, ranges)

# get hull area for all species
allhulls = hullarea.(allspecies)

# show the mean size of the climatic range in environmental space
overplot_pca_space(allcentroids, allhulls)

# get the geographical centroids for all species
gc = geocentroids.(allspecies)

# plot the location of all the centroids on the map
overplot_geo_space(gc, ranges; mask = sa_mask)

# divide the range sizes into quantiles for plotting
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

# Find the convex and concave hulls around all the points in pca space
## chull_all = LibGEOS.convexhull(points_to_geo(pca1, pca2))
bbox_all = (extrema(pca1), extrema(pca2))
cc = concave_hull(ConcaveHull.KDTree(vcat(pca1', pca2'), reorder = false), 40)
chull_all = GI.Polygon([GI.LinearRing([cc.vertices; [first(cc.vertices)]])])
# I am currently overwriting with a concave hull to see the effect
# This is the hull around all of the occupied space


# Now fit an ellipse to a species in pca space
plot_species_pca(rand(allspecies), 2)

# Plot 16 random species with occurrences in pca space and fitted ellipses
p = Plots.plot([
    plot_species_pca(rand(allspecies)) for i in 1:16]...
, size = (1200, 1200))
savefig(p, "16 species in pca space.png")

# Now control for the density of points in pca space by applying a 0.1 grid
weightmap = do_map(makeweights(pca1, pca2, 0.1), sa_mask)

p = Plots.plot([
    plot_species_pca(rand(allspecies), 1.5, weighted = true) for i in 1:16]...
, size = (1200, 1200))
savefig(p, "16 species controlling for point density.png")

# fit elliptical niches for all species
els = fitellipse.(allspecies)

# show patterns of ellipse area
ares = GeometryBasics.area.(els)
histogram(ares)
savefig("histogram of empirical ellipse areas.png")
Plots.scatter((el -> (el.center_x, el.center_y)).(els), marker_z = ares, ms = 3)
savefig("PCA centroids of empirical ellipses with area as color")

# repeat the plot with 
Plots.default(msw = 0, ms = 1, aspect_ratio = 1, seriescolor = cgrad(:Spectral, rev = true), legend = false, colorbar = true)
el_emp_point = [count(el -> in_ellipse(pt, el), els) for pt in zip(pca1, pca2)]
Plots.plot(
    Plots.scatter(pca1, pca2, marker_z = el_emp_point, title = "fitted ellipse overlap"), 
    Plots.scatter(pca1, pca2, marker_z = diversity[sa_mask], title = "empirical richness") 
)
savefig("empirical_ellipse_and_empirical_pca_richness.png")


# Create random ellipses with the empirical areas
ellipses = Ellipse[]
for harea in ares
    ## replace the line below with this to pick a random real grid cells as ellipse centers, rather than a random point
    # centerpoint = rand(1:length(pca1))
    # el = rand(Ellipse, pca1[centerpoint], pca2[centerpoint], area = area)

    el = rand(Ellipse, xlims = first(bbox_all), ylims = last(bbox_all), area = harea)
    ovrlp = harea != 0 ? overlap(el, chull_all) : 1
    failsafe = 0
    while ovrlp < 0.8 && (failsafe += 1) < 1000 
        el = rand(Ellipse, xlims = first(bbox_all), ylims = last(bbox_all), area = harea)
        ovrlp = harea != 0 ? overlap(el, chull_all) : 1
    end
    failsafe == 1000 && @show harea
    push!(ellipses, el)
end

# TODO possibly inside the loop above: grow a range based on the ellipse from a random point 
# Show 50 random ellipses
p = Plots.scatter(pca1, pca2, mc = :grey, ms = 1, msw = 0, aspect_ratio = 1, label = "")
for el in rand(ellipses, 50)
    Plots.plot!(p, el, label = "")
end
p
savefig(p, "50 random ellipses.png")

# plot the modelled and empirical richness
Plots.default(msw = 0, ms = 1, aspect_ratio = 1, seriescolor = cgrad(:Spectral, rev = true), legend = false, colorbar = true)
elpoint = [count(el -> in_ellipse(pt, el), ellipses) for pt in zip(pca1, pca2)]
Plots.plot(
    Plots.scatter(pca1, pca2, marker_z = elpoint, title = "ellipse overlap"), 
    Plots.scatter(pca1, pca2, marker_z = diversity[sa_mask], title = "empirical richness") 
)

savefig("modelled_ellipse_and_empirical_pca_richness.png")
