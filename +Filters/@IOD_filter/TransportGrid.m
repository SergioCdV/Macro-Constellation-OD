function [grid] = TransportGrid(obj, samples, mode, post_mode)
    % Preallocation 
    grid = samples;

    % Compute the rotation matrix from the previous mode to the new mode
    Q1 = QuaternionAlgebra.right_isoclinic(mode);
    Q2 = QuaternionAlgebra.right_isoclinic(post_mode);
    R = Q2 * Q1.'; 

    % Transport the grid 
    for i = 1:size(samples,2)
        grid(:,i) = R*samples(:,i);
    end
end