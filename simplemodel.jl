using Rasters, RasterDataSources
using Stencils
using MultivariateStats
using Distributions
using Plots
using Random
using StatsBase
using GLMakie
using Statistics

bioclim = RasterStack(CHELSA{BioClim}; lazy=true)
bioclim_sa = bioclim[X=-89 .. -33, Y=-57 .. 13]
bioclim_sa = aggregate(mean, replace_missing(bioclim_sa, NaN), 10)
sa_mask = boolmask(bioclim_sa.bio15)
Plots.plot(bioclim_sa.bio1)

# logged = RasterStack(Tuple(map(x -> x[1] > 11 ? log.(x[2]) : x[2], enumerate(values(replace_missing(bioclim_sa))))))
big_mat = reduce(vcat, map(enumerate(bioclim_sa)) do (i, A)
    permutedims(zscore(i > 11 ? log.(A[sa_mask] .+ 1) : A[sa_mask]))
end)
# histogram(big_mat'; bins=20, ticks=false, label=false, title=(1:19)')

model = fit(PCA, big_mat; maxoutdim=2)
pred = transform(model, big_mat)
pca1 = pred[1, :]
pca2 = pred[2, :]
p = Makie.scatter(collect(zip(pca1, pca2)); color=(:grey, 0.5))

map1 = fill(NaN, dims(bioclim_sa))
map2 = fill(NaN, dims(bioclim_sa))
map1[sa_mask] .= pca1
map2[sa_mask] .= pca2
Plots.plot(RasterStack(map1, map2))

# Rasterize

using Shapefile, GeoInterface, Extents#, ProgressMeter
const GI = GeoInterface
using GeometryOps
shapefiles = [
    "/home/raf/Data/Biodiversity/Distributions/Birds/batch_1.shp",
    "/home/raf/Data/Biodiversity/Distributions/Birds/batch_2.shp",
    "/home/raf/Data/Biodiversity/Distributions/Birds/batch_3.shp",
    "/home/raf/Data/Biodiversity/Distributions/Birds/batch_4.shp",
    "/home/raf/Data/Biodiversity/Distributions/Birds/batch_5.shp",
]
sa_geoms = map(shapefiles) do sf
    handle = Shapefile.Handle(shapefiles[2])
    filter(handle.shapes) do geom
        ext = GI.calc_extent(GI.trait(geom), geom)
        Extents.intersects(ext, Extents.extent(bioclim_sa))
    end
end |> Iterators.flatten |> collect

p = Makie.scatter(collect(zip(pca1, pca2)); color=(:grey, 0.1))
for i in 1:5
    species_mask = boolmask(sa_geoms[i]; to=sa_mask)
    Makie.scatter!(map1[species_mask], map2[species_mask]; 
        colormap=(:tab10, 0.1), color=i, colorrange=(1, 10)
    )
end
# save("climate_space_5_birds.png")

using GeometryBasics
simplified = GeometryOps.simplify(sa_geoms[1:100]; tol=0.0000001);
GeoInterface.convert(GeometryBasics, simplified)
rast = rasterize(count, sa_geoms; to=sa_mask, threaded=false)
Plots.plot(rast)


