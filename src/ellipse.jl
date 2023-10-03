using Random

# A basic Ellipse struct
struct Ellipse
    center_x::Float64
    center_y::Float64
    length::Float64
    width::Float64
    angle::Float64 # Given in radians
end

# check if the point `(xp, yp)`` is in an ellipse centered on (`x`,`y`) 
# with angle `an` (in radians), width `a` and height `b`
function in_ellipse((xp,yp), an, x, y, a, b)
    a = (cos(an) * (xp-x) + sin(an)*(yp-y))^2 / a^2
    b = (sin(an) * (xp-x) + cos(an)*(yp-y))^2 / b^2
    a+b < 1
end

function distance(point::Tuple, el::Ellipse)
    cosa = cos(el.angle)
    sina = sin(el.angle)
    rel_x = first(point) - el.center_x
    rel_y = last(point) - el.center_y
    a = (cosa * rel_x + sina * rel_y)^2 / el.length^2
    b = (sina * rel_x + cosa * rel_y)^2 / el.width^2
    a+b
end

in_ellipse(point, el::Ellipse) = distance(point, el) <= 1

rescale(x, min_x, max_x) = (max_x - min_x)*x + min_x

import Random.rand
function Random.rand(::Type{Ellipse}; 
    xlims=(0,1), ylims=(0,1), area = 1) 
    a = rand() * sqrt(pi)
    b = 1 / (pi * a)
    Ellipse(rescale(rand(), xlims...), rescale(rand(), ylims...), a * sqrt(area), b * sqrt(area), rand()π)
end