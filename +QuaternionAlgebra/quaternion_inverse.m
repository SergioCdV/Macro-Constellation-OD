

function [qi] = quaternion_inverse(q)
    qi = q;
    qi(1:3,1) = -qi(1:3,1);
end