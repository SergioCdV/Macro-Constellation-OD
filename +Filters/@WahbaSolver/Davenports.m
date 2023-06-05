
function  [q, Sigma] = Davenports(weights, b, r)
    % Attitude profile matrix
    B = b * diag(weights) * r.';

    % Davenports matrix 
    z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
    K = [B+B.'-trace(B) z; z.' trace(B)];

    % Optimal quaternion 
    [q, lambda] = eigs(K,1);
    A = QuaternionAlgebra.Quat2Matrix(q);

    % Compute the scaling parameter 
    c = sum(weights);

    % Covariance matrix
    Sigma = c * ((lambda * eye(3) - B * A) \ eye(3));
end