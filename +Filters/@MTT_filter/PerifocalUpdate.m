
function [plane, Sigma] = PerifocalUpdate(obj, weights, particles)
    % Reliability weights
    c = sum(weights,2);
    weights = weights / c;
    
    % Final mean estimation 
    [q, Sigma_q] = QuaternionAlgebra.AverageQuat(ones(1,size(particles,2)), particles(1:4,:));

    % Action estimation 
    beta = sum(weights .* particles(5:end,:),2);

    % Covariance estimation
    Sigma_b = (weights.*particles(5:7,:)-beta*weights) * (weights.*particles(5:7,:)-beta*weights).'/(1-sum(weights.^2,2));
    Sigma = blkdiag(Sigma_q, Sigma_b);
    Sigma = 0.5 * (Sigma + Sigma.') + obj.PD_tol * eye(size(Sigma));

    % Uscented transform 
%     UF = Filters.UKF('UKF-A', 2, 1e-1, 0);
%     UF.StateDim = 6;
%     UF = UF.Init();
%     sigma = UF.sigma_points(x, Sigma);
%     X = sum(UF.W(1,:).*sigma,2);

    X = zeros(6,1);

    % Mean update 
    plane = [q; beta;];                                                  % Mean state
    plane(1:4) = QuaternionAlgebra.exp_map([X(1:3); 0], plane(1:4));     % Perifocal attitude update
    plane(5:end) = plane(5:end) + X(4:end);                              % Action set update

    % Covariance update
%     Sigma_c = (sigma-X)*diag(UF.W(2,:))*(sigma-X).';
%     Sigma = Sigma_c;

    % Numerical conditioning
%     Sigma = Sigma / (length(weights)-1);
%     Sigma = 0.5 * (Sigma + Sigma.') + obj.PD_tol * eye(size(Sigma));
end