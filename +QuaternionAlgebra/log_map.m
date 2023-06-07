

function [l] = log_map(x, v)
    cos_alpha = v.' * x;

    if (cos_alpha ~= 1)
        l = (x - cos_alpha * v)  * acos(cos_alpha) / sqrt(1-cos_alpha^2);
    else
        l = [0;0;0;1];
    end

    l = [x(1:3)/norm(x(1:3)) * atan2(norm(x(1:3)), x(4)); 0];
end