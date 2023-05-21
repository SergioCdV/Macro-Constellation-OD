

function [e] = exp_map(x, v)
    n = norm(x(1:3)); 

    if (n)
        e = x / n * sin(n) + cos(n) * v;
    else
        e = v;
    end
end