
function  [q, Sigma] = ESOQ(obj, weights, b, r)
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
    H = B - lambda * eye(3);

    q = zeros(5, size(H,2));
    for i = 1:size(H,2)
        index = [1:i-1 i+1:size(H,2)];
        h = H(index,index);
        q(1,i) = det(h);
        q(1+i,i) = -det(h);
        q(2:i-1,i+1:end) = adj(h) * H([1:i-1 i+1:size(H,2)],i);
    end

    [~, index] = sort(q(1,:)); 
    q = q(:,index(end));
    q = q / norm(q);
    A = QuaternionAlgebra.Quat2Matrix(q);

    % Compute the scaling parameter 
    c = sum(weights);

    % Covariance matrix
    Sigma = c * ((lambda * eye(3) - B * A) \ eye(3));
end

%% Auxiliary functions 
