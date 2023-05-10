
function [Q] = left_isoclinic(q)
    Q = [q(4) * eye(3)+QuaternionAlgebra.hat_map(q(1:3)) q(1:3); -q(1:3).' q(4)];
end