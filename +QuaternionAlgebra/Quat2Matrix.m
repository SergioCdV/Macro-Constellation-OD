

function [A] = Quat2Matrix(q)
    S = QuaternionAlgebra.hat_map(q(1:3,1));
    A = (q(4,1)^2 - q(1:3,1).' * q(1:3,1)) * eye + 2 * q(1:3,1) * q(1:3,1).' - 2 * q(4) * S;
end