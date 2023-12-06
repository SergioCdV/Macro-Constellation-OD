
function [v] = RotateVector(q, v)
    aux = [v; 0]; 
    inv_quat = QuaternionAlgebra.quaternion_inverse(q);
    aux = QuaternionAlgebra.right_isoclinic(aux) * inv_quat;
    aux = QuaternionAlgebra.right_isoclinic(q) * aux;
    v = aux(1:3,1);
end