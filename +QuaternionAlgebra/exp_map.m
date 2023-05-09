

function [e] = exp_map(x, v)
    e = cos(norm(x)) * v + x * sin(norm(x)) / norm(x);
end