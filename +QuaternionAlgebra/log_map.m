

function [l] = log_map(x, v)
    cos_alpha = v.' * x;
    l = (x - cos_alpha * v)  * acos(cos_alpha) / sqrt(1-cos_alpha^2);
end