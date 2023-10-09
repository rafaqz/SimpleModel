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
    on_domain(pt, domain) || error("start point not on domain")
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

function random_point_on_ellipse(el::Ellipse, x = pca1, y = pca2; maxiter = 1e6)
    iter = 0
    while (iter += 1) < maxiter
        pt = rand(1:length(x))
        in_ellipse((x[pt], y[pt]), el) && return pt
    end
    error("Did not find a point on the ellipse in $maxiter tries")
end

map_ellipse(el::Ellipse, pca1=pca1, pca2=pca2, mask = sa_mask) = do_map([in_ellipse(pt, el) for pt in zip(pca1, pca2)], mask; missingval = false)

const sa_inds = collect(Iterators.product(1:size(sa_mask, 1), 1:size(sa_mask, 2)))[sa_mask]

function make_continuous_range(el, pca1=pca1, pca2=pca2, mask=sa_mask)
    domain = map_ellipse(el)
    i = random_point_on_ellipse(el, pca1, pca2)
    pt = sa_inds[i]
    @show i, pt
    growrange(pt, domain)
end

# TODO weird bug here
newdiv = reduce(+, make_continuous_range(el) for el in els)

myel = rand(els)
Plots.plot(
    Plots.plot(map_ellipse(myel)),
    Plots.plot(make_continuous_range(myel))
)