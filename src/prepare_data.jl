using Rasters
using RasterDataSources
using MultivariateStats
using DataFrames
using ArchGDAL
using GeoInterface; const GI = GeoInterface
using LibGEOS
using Stencils
using Extents
using Shapefile

# Let's load the environment data

function prepare_environment(datadir)
    ENV["RASTERDATASOURCES_PATH"] = joinpath(datadir, "Rasterdatasources")
    # datadir = "/home/raf/Data/Biodiversity/Distributions"

    # download the bioclim variables for south america

    # bioclim = RasterStack(CHELSA{BioClim}; lazy=true)
    # bioclim_sa = bioclim[X=-89 .. -33, Y=-57 .. 13]
    # bioclim_sa = Rasters.aggregate(mean, replace_missing(bioclim_sa, NaN), 10)

    # Use WorldClim for now for speed
    bioclim = RasterStack(WorldClim{BioClim}; lazy=true)
    bioclim[X=-89 .. -33, Y=-57 .. 13]
    # bioclim_sa = Rasters.aggregate(mean, replace_missing(bioclim_sa, NaN), 10)
end

# Fit a PCA model to the climate and extract the two primary components
function do_pca(bioclim_sa, sa_mask)
    # logged = RasterStack(Tuple(map(x -> x[1] > 11 ? log.(x[2]) : x[2], enumerate(values(replace_missing(bioclim_sa))))))
    big_mat = reduce(vcat, map(enumerate(bioclim_sa)) do (i, A)
        permutedims(zscore(i > 11 ? log.(A[sa_mask] .+ 1) : A[sa_mask]))
    end)
    # histogram(big_mat'; bins=20, ticks=false, label=false, title=(1:19)')

    model = fit(PCA, big_mat; maxoutdim=2)
    pred = MultivariateStats.transform(model, big_mat)

    pred[1, :], pred[2, :], model
end

# Load the bird shapefiles and pick the ones in South America

function loadranges(data::String, batches::Int, mask, datadir)
    shapefiles = [joinpath(datadir, data, "batch_$i.shp") for i in 1:batches]
    reduce(vcat, map(shapefiles) do sf
        df = DataFrame(Shapefile.Table(sf))
        filter(df) do row
            ext = GI.calc_extent(GI.trait(row.geometry), row.geometry)
            Extents.intersects(ext, Extents.extent(mask))
        end
    end)
end

