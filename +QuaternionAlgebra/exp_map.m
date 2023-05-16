

function [e] = exp_map(x, v)
    n = norm(x); 

    if (n)
        e = x / n * sin(n) + cos(n) * v;
    else
        e = [0;0;0;1];
    end
end