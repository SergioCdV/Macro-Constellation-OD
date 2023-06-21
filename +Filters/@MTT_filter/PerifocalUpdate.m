
function [plane, Sigma] = PerifocalUpdate(obj, weights, particles)
    % Reliability weights
    c = sum(weights,2);
    if (c)
        weights = weights / c;
    end
    
    % Final mean estimation 
    [q, ~] = QuaternionAlgebra.AverageQuat(ones(1,size(particles,2)), particles(1:4,:));

    a = zeros(3,size(particles,2));
    for i = 1:size(particles,2)
        dq = QuaternionAlgebra.left_isoclinic( particles(1:4,i) ).' * q;
        a(:,i) = QuaternionAlgebra.MPR2Quat(1, 1, dq, false);
    end

    % Action estimation 
    am = sum(weights .* a, 2);
    beta = sum(weights .* particles(5:end,:),2);

    % Covariance estimation
    x = [a; particles(5:7,:)];
    mu = [am; beta];
    Sigma = (weights.*x-mu*weights) * (weights.*x-mu*weights).'/(1-sum(weights.^2,2));
    Sigma = 0.5 * (Sigma + Sigma.') + obj.PD_tol * eye(size(Sigma,1));

    % Final output 
    plane = [q; beta];
end