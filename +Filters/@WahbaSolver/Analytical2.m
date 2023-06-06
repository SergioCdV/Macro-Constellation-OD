
function [q, Sigma] = Analytical2(obj, weights, b, r)
    % Compute the unnormalized quaternion
    q(1:3,1) = cross(b(:,1)-r(:,1), b(:,2)-r(:,2));
    q(4,1) = dot(b(:,2),r(:,1))- dot(b(:,1),r(:,2));

    % Normalization 
    q = q / norm(q);

    % Covariance
    Sigma = 1e-10 * eye(3);
end