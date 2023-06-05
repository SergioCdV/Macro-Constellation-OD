

function [CorrPlane, CorrPlaneSigma] = PlaneCorrection(obj, PlaneEstimator, planes, Sigma, sigma_points, Posterior)
    % Box constraints 
    box(1,:) = [1 1.5];
    box(2,:) = [sqrt(1-0.01^2) 1];
    box(3,:) = [-1 1];

    % Covariance estimation
    omega = Posterior(5:7,:);

    Cov = zeros(3 * size(omega,2));
    for i = 1:size(Posterior,2)
        Sigma_o = reshape(Posterior(9:end,i), [7 7]);
        Cov(1+3*(i-1):3*i, 1+3*(i-1):3*i) = Sigma_o(4:6,4:6);
    end

    % UKF step     
    PlaneEstimator.R = Cov;
    PlaneEstimator = PlaneEstimator.AssignObservationProcess(3 * size(omega,2), @(state)MeasurementModel(size(omega,2), state));
    omega = reshape(omega, [], 1);
    [CorrPlane, CorrPlaneSigma, ~, ~] = PlaneEstimator.CorrectionStep(sigma_points, planes, Sigma, omega);

    % ADMM step 
    X = CorrPlane(4:6,1);
    for k = 1:size(X,2)
        % Projection of the Delaunay action 
        if (X(1,k) < box(1,1))
            X(1,k) = box(1,1);
        elseif (X(1,k) > box(1,2))
            X(1,k) = box(1,2);
        end
        
        % Projection of the angular momentum 
        if (X(2,k) / X(1,k) < box(2,1))
            X(2,k) = X(1,k) * box(2,1);
        elseif (X(2,k) / X(1,k) > box(2,2))
            X(2,k) = X(1,k) * box(2,2);
        end
    
        % Projection of the nodal angular momentum 
        if (X(3,k) / X(2,k) < box(3,1))
            X(3,k) = X(2,k) * box(3,1);
        elseif (X(3,k) / X(2,k) > box(3,2))
            X(3,k) = X(2,k) * box(3,2);
        end
    end
    
    % Final output 
    CorrPlane = [CorrPlane(7:10,1); X];
end

%% Auxiliary function 
function [y] = MeasurementModel(dim, state)
    % Preallocation 
    y = zeros(3 * dim, size(state,2));

    % Measurement of the angulat velocity
    for i = 1:size(state,2)
        y(:,i) = repmat(state(5:7,1), dim, 1);
    end
end