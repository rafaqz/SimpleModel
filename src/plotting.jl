using Plots
using GLMakie
using MultivariateStats
using GeoInterfaceMakie
using GeoInterfaceRecipes


# visualize a PCA biplot

function biplot(pca1, pca2, model, labels = [""])
    p = Makie.scatter(collect(zip(pca1, pca2)); color=(:grey, 0.5))
    for i in 1:19
        Makie.lines!([(0,0), 15 .* (projection(model)[i,:]...,)], color = :red)
        Makie.text!((15 .* projection(model)[i,:])..., text = labels[i])
    end
    p
end

# This plots the distribution of a random species on the map
function plot_distribution(species, mask = sa_mask)
    Plots.plot(maske, fill = :grey)
    Plots.plot!(get_speciesmask(rand(allspecies)))
end

# Creates two side-by-side plots, one in geographic space, the other in climate space
# Shows the occurrence points of all the species defined by speciesnames
function plot_species(speciesnames; mask = sa_mask, pca1 = pca1, pca2 = pca2, pca_maps = pca_maps)
    f = Figure(resolution = (1500, 700))
    a = Axis(f[1,1], aspect = DataAspect())
    b = Axis(f[1,2], aspect = DataAspect())
    Makie.plot!(a, mask, colormap = "Greys")
    Makie.scatter!(b, collect(zip(pca1, pca2)); markersize = 0.1, color=(:grey, 0.5))
    colsize!(f.layout, 1, Auto(0.73))
    rich = fill(NaN, dims(sa_mask))

    for (i, spec) in enumerate(speciesnames)
        species_mask = get_speciesmask(spec)
        rich[species_mask] .= i
        Makie.scatter!(b, pca_maps.pca1[species_mask], pca_maps.pca2[species_mask]; 
            markersize = 0.2, label = spec,
            colormap=(:tab10, 0.5), color=i, colorrange=(1, 10)
        )
    end

    Makie.plot!(a, rich, colormap=(:tab10, 0.5), colorrange=(1, 10))
  #  leg = Makie.Legend(f[1,3], a, framevisible = false, colorrange = (1,10)) # TODO: there was a problem with the legend code
  #  for i in 1:length(speciesnames)
  #      leg.entrygroups[][1][2][i].elements[1].markersize = 20
  #  end
    f
end


# Overplot some points on the pca background
function overplot_pca_space(points, color; pca1=pca1, pca2=pca2)
    p = Figure()
    a = Axis(p[1,1], aspect = DataAspect())
    Makie.scatter!(a, collect(zip(pca1, pca2)); markersize = 0.3, color=(:grey, 0.5))
    Makie.scatter!(a, points; markersize = 4, color = color, label = "rangesizes", colormap = Reverse(:Spectral))
    Makie.Colorbar(p[1,2], colormap = Reverse(:Spectral))
    p
end

# Overplot some points on the map of south america
function overplot_geo_space(points, col; mask = sa_mask)
    f = Figure()
    a = Axis(f[1,1], aspect = DataAspect())
    Makie.plot!(a, mask, colormap = :Greys)
    Makie.scatter!(points, color = col)
    f
end


