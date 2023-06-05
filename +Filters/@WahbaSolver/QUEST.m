
function  [q, Sigma] = QUEST(weights, b, r)
    % Attitude profile matrix
    B = b * diag(weights) * r.';

    % Davenports matrix 
    z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
    H = B + B.';

    delta = det(H);
    sigma = trace(B);
    k = trace(adj(H));
    a = sigma^2-k;
    b = sigma^2+z.'*z;
    c = delta + z.' * H * z;
    d = z.' * (H * H) * z;

    % Optimal quaternion 
    fun = @(lambda)(lambda^4-(a+b)*lambda^2-c*lambda+(a*b+c*sigma-d));
    lambda = fsolve(fun, 1);
    p = ((lambda + sigma) * eye(3) -H) \ z;
    q = [p; 1]; 
    q = q / norm(q);
    A = QuaternionAlgebra.Quat2Matrix(q);

    % Compute the scaling parameter 
    c = sum(weights);

    % Covariance matrix
    Sigma = c * ((lambda * eye(3) - B * A) \ eye(3));
end

%% Auxiliary functions 
