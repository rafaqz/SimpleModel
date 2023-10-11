using Plots
using GLMakie
using Rasters
using JLSO

include("objects.jl")
include("plotting.jl")
include("ellipse.jl")
include("simplemodel.jl")


###--- First we load all the data into two objects. This takes a while if first time, 
# so we use JLSO to cache

datadir = "/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data"
obj = try
    JLSO.load(joinpath(datadir, "processed_objects.jls"))
catch
    include("prepare_data.jl")
    obj = prepare_data(datadir)
end
spec = obj[:spec] # see the "objects.jl" file for an explanation of these two objects
env = obj[:env]

###--- Exploratory data analysis

# Count total diversity
diversity = reduce(+, spec.ranges)
Plots.heatmap(diversity, color = cgrad(:Spectral, rev = true))
savefig("figures/empirical_richness.png")


# Plot the diversity in climate space
f = Figure()
a, s = Makie.scatter(f[1,1], collect(zip(env.pca1, env.pca2)); markersize = 0.1, color = diversity[env.mask], colormap = cgrad(:Spectral, rev = true))
Colorbar(f[1,2],s)
display(f)


# plot some random species and look at their distribution
randspecs = rand(spec.names, 5)
plot_species(randspecs, spec, env)


## how are ranges shaped?
cors, xrange, yrange = find_range_shapes(spec, env)
histogram(cors)
Plots.scatter(xrange, yrange)


###--- Patterns of range and niche size in climate space

# get environmental centroids and range sizes for all species
allcentroids = vec(get_centroid.(spec.ranges, (env,)))
rangesizes = vec(count.(spec.ranges))

# This figure shows the geographic range size of ranges on the environmental centroids (in environment spaecs)
overplot_pca_space(allcentroids, rangesizes, env)

# get hull area for all species
allhulls = [hullarea(sp, spec, env) for sp in spec.names]

# show the mean size of the climatic range in environmental space
overplot_pca_space(allcentroids, allhulls, env)

# get the geographical centroids for all species
gc = geocentroids.(spec.names, (spec,))

# plot the location of all the centroids on the map
overplot_geo_space(gc, rangesizes, env)

# divide the range sizes into quantiles for plotting
rangequants = asquantile(rangesizes, 4)

# and plot them
f = Figure()
aspect = DataAspect()
axs = [Axis(f[1,1]; aspect), Axis(f[2,1]; aspect), Axis(f[1,2]; aspect), Axis(f[2,2]; aspect)]
for i in 1:4
    inds = findall(==(i), rangequants)
    Makie.plot!(axs[i], env.mask, colormap = :Greys)
    Makie.scatter!(axs[i], gc[inds], markersize = 2)
end
f

# Now fit an ellipse to a species in pca space
plot_species_pca(rand(spec.names), spec, env, 2)

# Plot 16 random species with occurrences in pca space and fitted ellipses
p = Plots.plot([
    plot_species_pca(rand(spec.names), spec, env, 2) for i in 1:16]...
, size = (1200, 1200))
savefig(p, "figures/16 species in pca space.png")

# Now control for the density of points in pca space by applying a 0.1 grid
# First map the options, e.g. binsizes from 0.1 to 0.4
map_binsize(binsize) = Plots.heatmap(1 ./ do_map(makeweights(env.pca1, env.pca2, binsize), env.mask), color = cgrad(:Spectral, rev = true), title = "binsize = $binsize")
Plots.plot([map_binsize(bs) for bs in 0.1:0.1:0.4]..., size = (800, 800))
savefig("figures/binsizes.png")

# we conclude that we need one at 0.2
const weightmap = do_map(makeweights(env.pca1, env.pca2, 0.2), env.mask)

p = Plots.plot([
    plot_species_pca(rand(spec.names), spec, env, 1.5; weightmap) for i in 1:16]...
, size = (1200, 1200))
savefig(p, "figures/16 species controlling for point density.png")

# fit elliptical niches for all species
emp_ellipses = [fitellipse(name, spec, env, 2) for name in spec.names]

# show patterns of ellipse area
ares = GeometryBasics.area.(emp_ellipses)
histogram(ares)
savefig("figures/histogram of empirical ellipse areas.png")
Plots.scatter((el -> (el.center_x, el.center_y)).(emp_ellipses), marker_z = ares, ms = 3)
savefig("figures/PCA centroids of empirical ellipses with area as color")


# repeat the plot with the 
Plots.default(msw = 0, ms = 1, aspect_ratio = 1, seriescolor = cgrad(:Spectral, rev = true), legend = false, colorbar = true)
el_emp_point = [count(el -> in_ellipse(pt, el), emp_ellipses) for pt in zip(env.pca1, env.pca2)]
Plots.plot(
    Plots.scatter(env.pca1, env.pca2, marker_z = el_emp_point, title = "fitted ellipse overlap"), 
    Plots.scatter(env.pca1, env.pca2, marker_z = diversity[env.mask], title = "empirical richness") 
)
savefig("figures/empirical_ellipse_and_empirical_pca_richness.png")

  
# Create random ellipses with the empirical areas
rand_ellipses = [sample_ellipse(harea, env; on_real_point = true) for harea in ares]

# Show 50 random ellipses
p = Plots.scatter(env.pca1, env.pca2, mc = :grey, ms = 1, msw = 0, aspect_ratio = 1, label = "")
for el in rand(rand_ellipses, 50)
    Plots.plot!(p, el, label = "")
end
p
savefig(p, "figures/50 random ellipses.png")

# plot the modelled and empirical richness
Plots.default(msw = 0, ms = 1, aspect_ratio = 1, seriescolor = cgrad(:Spectral, rev = true), legend = false, colorbar = true)
elpoint = [count(el -> in_ellipse(pt, el), rand_ellipses) for pt in zip(env.pca1, env.pca2)]
Plots.plot(
    Plots.scatter(env.pca1, env.pca2, marker_z = elpoint, title = "ellipse overlap"), 
    Plots.scatter(env.pca1, env.pca2, marker_z = diversity[env.mask], title = "empirical richness") 
)

savefig("figures/modelled_ellipse_and_empirical_pca_richness.png")

Plots.heatmap(do_map(el_emp_point, env.mask), color = cgrad(:Spectral, rev = true))
savefig("figures/richness_on_empirical_ellipses.png")

Plots.heatmap(do_map(elpoint, env.mask), color = cgrad(:Spectral, rev = true))
savefig("figures/richness_on_random_ellipses.png")

model_ranges = RasterSeries([make_continuous_range(el, env) for el in emp_ellipses], (; name = spec.names))
newdiv = reduce(+, model_ranges)
Plots.heatmap(newdiv, color = cgrad(:Spectral, rev = true))
savefig("figures/richness from patches in empirical ellipses.png")

# Find and plot richnes at the 1 degree lat/long scale
model_coarse = reduce(+, Rasters.aggregate.(any, model_ranges, 6))
emp_coarse = reduce(+, Rasters.aggregate.(any, spec.ranges, 6))

Plots.default(fillcolor = cgrad(:Spectral, rev = true))
Plots.plot(
    Plots.heatmap(model_coarse, title = "modelled 1 degree richness"),
    Plots.heatmap(emp_coarse, title = "empirical 1 degree richness"), size = (900, 500)
)
savefig("figures/coarse richness.png")
