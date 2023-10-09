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

# Find the right rotation for the PCA - copied from factorloadingmatrices (varimax)
function vmax(A::AbstractMatrix{TA}; gamma = 1.0, minit = 20, maxit = 1000,
    reltol = 1e-12) where TA
    d, m = size(A)
    m == 1 && return A
    # Warm up step: start with a good initial orthogonal matrix T by SVD and QR
    T = Matrix{TA}(I, m, m)
    B = A * T
    L,_,M = svd(A' * (d*B.^3 - gamma*B * Diagonal(sum(B.^2, dims = 1)[:])))
    T = L * M'
    if norm(T - Matrix{TA}(I, m, m)) < reltol
        T,_ = qr(randn(m,m)).QT
        B = A * T
    end

    # Iteration step: get better T to maximize the objective (as described in Factor Analysis book)
    D = 0
    for k in 1:maxit
        Dold = D
        L,s,M = svd(A' * (d*B.^3 - gamma*B * Diagonal(sum(B.^2, dims = 1)[:])))
        T = L * M'
        D = sum(s)
        B = A * T
        if (abs(D - Dold)/D < reltol) && k >= minit
            break
        end
    end
    return T
end



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
    pred2 = -(pred' * vmax(loadings(model))) # the minus here and below are just a transformation to have high values top right

    pred2[:, 1], pred2[:, 2], -(loadings(model) * vmax(loadings(model)))
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

