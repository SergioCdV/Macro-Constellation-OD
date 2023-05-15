
function [Y] = quadrature(obj, x)
    Y = dot(obj.dt(2,:) .* obj.W, x, 2);
end