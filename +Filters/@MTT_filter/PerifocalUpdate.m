
function [plane, Sigma] = PerifocalUpdate(obj, weights, particles)
    % Reliability weights
    c = sum(weights,2);
    weights = weights / c;
    
    % Final mean estimation 
    [q, Sigma_q] = QuaternionAlgebra.AverageQuat(ones(1,size(particles,2)), particles(1:4,:));

    a = zeros(3,size(particles,2));
%     for i = 1:size(particles,2)
%         dq = QuaternionAlgebra.left_isoclinic( particles(1:4,i) ) * q;
%         aux = QuaternionAlgebra.log_map(dq, [0;0;0;1]);
%         a(:,i) = aux(1:3);
%     end

    % Action estimation 
    beta = sum(weights .* particles(5:end,:),2);

    % Covariance estimation
    x = [a; particles(5:7,:)];
    mu = [zeros(3,1); beta];
    Sigma = (weights.*x-mu*weights) * (weights.*x-mu*weights).'/(1-sum(weights.^2,2));
    Sigma(1:3,1:3) = Sigma_q;
    Sigma = 0.5 * (Sigma + Sigma.') + obj.PD_tol * eye(size(Sigma,1));

    % Final output 
    plane = [q; beta];
end