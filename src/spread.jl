const nbh = Tuple((x,y) for x in -1:1, y in -1:1 if !(x == y == 0))

on_domain(pt, domain) = min(pt...) > 0 && first(pt) <= size(domain, 1) && last(pt) <= size(domain, 2) && domain[pt...]

function add_point!(georange, potentials, pt, domain)
    georange[pt...] = true
    for nb in nbh
        newpt = nb .+ pt
        if on_domain(newpt, domain) && !georange[newpt...]
            push!(potentials, newpt)
        end
    end
end

function growrange(start, domain)
    on_domain(start, domain) || error("start point not on domain")
    georange = fill(false, dims(domain), missingval = false)
    potentials = Set([start])
    failsafe = 0; max_iter = prod(size(domain))
    while !isempty(potentials)
        (failsafe += 1) > max_iter && error("stuck")
        pt = pop!(potentials)
        add_point!(georange, potentials, pt, domain)
    end
    georange
end

function random_point_on_ellipse(el::Ellipse, x, y; maxiter = 1e6)
    iter = 0
    while (iter += 1) < maxiter
        pt = rand(1:length(x))
        in_ellipse((x[pt], y[pt]), el) && return pt
    end
    error("Did not find a point on the ellipse in $maxiter tries")
end
random_point_on_ellipse(el::Ellipse, env::Environment; maxiter = 1e6) = random_point_on_ellipse(el, env.pca1, env.pca2; maxiter)

map_ellipse(el::Ellipse, env::Environment) = do_map([in_ellipse(pt, el) for pt in zip(env.pca1, env.pca2)], env)

function make_continuous_range(el, env::Environment)
    domain = map_ellipse(el, env)
    GeometryBasics.area(el) == 0 && return domain
    i = random_point_on_ellipse(el, env.pca1, env.pca2)
    pt = env.inds[i]
    growrange(pt, domain)
end



model_ranges = RasterSeries([make_continuous_range(el, env) for el in emp_ellipses], (; name = spec.names))
newdiv = reduce(+, model_ranges)
Plots.heatmap(newdiv, color = cgrad(:Spectral, rev = true))
savefig("figures/richness from patches in empirical ellipses.png")


model_coarse = reduce(+, Rasters.aggregate.(any, model_ranges, 6))
emp_coarse = reduce(+, Rasters.aggregate.(any, spec.ranges, 6))

Plots.default(fillcolor = cgrad(:Spectral, rev = true))
Plots.plot(
    Plots.heatmap(model_coarse, title = "modelled 1 degree richness"),
    Plots.heatmap(emp_coarse, title = "empirical 1 degree richness"), size = (900, 500)
)
savefig("figures/coarse richness.png")



myel = rand(els)
Plots.plot(
    Plots.plot(map_ellipse(mye, env)),
    Plots.plot(make_continuous_range(myel, env))
)
