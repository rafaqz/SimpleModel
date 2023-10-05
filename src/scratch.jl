
# Test the ellipse code #TODO move to plotting and script
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

as, bs, es = Float64[], Float64[], Float64[]
for i in 1:100
    e = rand(Ellipse)
    a, b = e.length, e.width
    push!(es, a/b)
    push!(as, a)
    push!(bs, b)
end
Plots.scatter(as, bs)
Plots.histogram(es)



# plot a species on the density of points in climate space
background = fit(Histogram, (pca1, pca2),  (-10:0.3:10, -10:0.3:10))
spec = rand(allspecies)
x,y = get_climate(allranges[At(spec)])
sp = fit(Histogram, (x,y),  (-10:0.3:10, -10:0.3:10))
Plots.plot(Plots.plot(sp, aspect_ratio = 1), Plots.plot(background, aspect_ratio = 1), layout = 2)

# relative occupancy of climate space
rel = sp.weights ./ background.weights
replace!(rel, 0 => NaN)
Plots.heatmap(rel')




# not really necessary
#using GeometryBasics
#simplified = GeometryOps.simplify(sa_geoms[1:100]; tol=0.0000001);
#GeoInterface.convert(GeometryBasics, simplified)
#rast = rasterize(count, sa_geoms; to=sa_mask, threaded=false)
#Plots.plot(rast)


el = rand(ellipses)
els = [in_ellipse(pt, el) for pt in zip(pca1, pca2)]
elmap = do_map(els, sa_mask; missingval = false)


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
