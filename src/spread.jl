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

# TODO weird bug here
newdiv = reduce(+, make_continuous_range(el, env) for el in emp_ellipses)

myel = rand(els)
Plots.plot(
    Plots.plot(map_ellipse(mye, env)),
    Plots.plot(make_continuous_range(myel, env))
)
for i in 1:length(emp_ellipses)
    @show i
    make_continuous_range(emp_ellipses[i], env)
end