

function [l] = log_map(x, v)
    alpha = acos(v.' * x);
    l = (x - cos(alpha) * v)  * alpha / sin(alpha);
end