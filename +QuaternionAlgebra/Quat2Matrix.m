

function [A] = Quat2Matrix(q)
    A = (q(4,1)^2 - dot(q(1:3,1), q(1:3,1))) * eye(3) + 2 * q(1:3,1) * q(1:3,1).' - 2 * q(4,1) * QuaternionAlgebra.hat_map(q(1:3,1));
end