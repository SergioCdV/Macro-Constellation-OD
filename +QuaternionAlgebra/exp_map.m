

function [e] = exp_map(x, v)
    n = norm(x(1:3)); 

    if (n)
        e = [x(1:3) / n * sin(n); cos(n)];
    else
        e = [0;0;0;1];
    end
end