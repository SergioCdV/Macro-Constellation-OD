
function [q, Sigma] = AverageQuat(w, dq)

    % Compute the average quaternion through the Markeley 2007 method
    Q = w .* dq; 
    M = Q * Q.';

    % Get the eigenvector corresponding to the largest eigenvalue 
    [q,~] = eigs(M,1);

    % Compute the covariance matrix
    R = zeros(4);

    Q = QuaternionAlgebra.left_isoclinic(q);
    for i = 1:size(dq,2)
        Qi = QuaternionAlgebra.left_isoclinic(dq(:,i));
        R = R + w(i) * Qi(:,1:3) * Qi(:,1:3).';
    end
    Sigma = (Q(:,1:3).' * R * Q(:,1:3))^(-1);
end