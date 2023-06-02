

function [A] = Quat2Matrix(q)
    A = QuaternionAlgebra.left_isoclinic(q).' * QuaternionAlgebra.right_isoclinic(q);
end