using Plots
using GLMakie
using Rasters

include("prepare_data.jl")
include("plotting.jl")
include("ellipse.jl")
include("simplemodel.jl")

###--- First we load all the data

## Get the environmental data
bioclim_sa = prepare_environment("/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data")
sa_mask = boolmask(bioclim_sa.bio15)

# and visualize
Plots.plot(bioclim_sa.bio1)

## get the PCA
pca1, pca2, model = do_pca(bioclim_sa, sa_mask)

# and convert the results to raster
pca_maps = RasterStack((pca1=do_map(pca1, sa_mask), pca2=do_map(pca2, sa_mask)))

# and visualize
Plots.plot(pca_maps)
biplot(pca1, pca2, model, string.(names(bioclim_sa)))

# Get the bird data (this takes time)
sa_geoms = loadranges("Birds", 5, sa_mask, "/Users/cvg147/Library/CloudStorage/Dropbox/Arbejde/Data")
# names of all species
allspecies = unique(sa_geoms.sci_name)

###--- Exploratory data analysis

# Count total diversity
diversity = rasterize(count, sa_geoms; to=sa_mask, boundary=:touches)
diversity .*= sa_mask
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
histogram(cors)
Plots.scatter(xrange ./ yrange)


###--- Patterns of range and niche size in climate space

# get environmental centroids and range sizes for all species
allcentroids = get_centroid.(get_speciesmask.(allspecies))
ranges = count.(get_speciesmask.(allspecies))

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

# and plot them #TODO simplify
f = Figure()
aspect = DataAspect()
axs = [Axis(f[1,1]; aspect), Axis(f[2,1]; aspect), Axis(f[1,2]; aspect), Axis(f[2,2]; aspect)]
for i in 1:4
    inds = findall(==(i), rangequants)
    Makie.plot!(axs[i], sa_mask, colormap = :Greys)
    Makie.scatter!(axs[i], gc[inds], markersize = 2)
end
f

chull_all = LibGEOS.convexhull(points_to_geo(pca1, pca2))
bbox_all = (extrema(pca1), extrema(pca2))

ellipses = Ellipse[]
for area in allhulls
    el = rand(Ellipse, xlims = first(bbox_all), ylims = last(bbox_all), area = area)
    ovrlp = area != 0 ? overlap(el, chull_all) : 1
    while ovrlp < 0.9
        el = rand(Ellipse, xlims = first(bbox_all), ylims = last(bbox_all), area = area)
        ovrlp = area != 0 ? overlap(el, chull_all) : 1
    end
    push!(ellipses, el)
end

# TODO possibly inside the loop above: grow a range based on the ellipse from a random point 

p = Plots.scatter(pca1, pca2, mc = :grey, ms = 1, msw = 0, aspect_ratio = 1, label = "")
for el in rand(ellipses, 50)
    Plots.plot!(p, el, label = "")
end
p