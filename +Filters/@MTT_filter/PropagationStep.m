
function [PropPrior, sigma_points] = PropagationStep(obj, last_epoch, new_epoch, Estimator, Prior)
    % Constants 
    step = new_epoch - last_epoch;      % Time step
    step = 86400 / obj.Tc * step;       % Dimensionalizing
    pos = 8; 

    % Preallocation 
    particles = Prior(1:end-1,:);       % Prior states
    weights = Prior(end,:);             % Prior weights
    
    % Reduce the dimensionality dynamics 
    REstimator = Estimator;
    REstimator.StateDim = 4; 
    REstimator.Q = REstimator.Q(pos-4:pos-1, pos-4:pos-1);

    sigma_points = zeros(pos * (2*REstimator.StateDim+1), size(particles,2));

    % Perform the propagation through the KF step
    REstimator = REstimator.Init().AssignStateProcess(REstimator.StateDim, @(M, step)dynamics(obj.epsilon, M, step));
    for i = 1:size(particles,2)
        % Conditioning
        Sigma = reshape(particles(pos+1:end,i), [pos-1 pos-1]); 
        mu = particles(pos-3:pos,i);
        sigma = Sigma(4:end,4:end); 

        [~, flag] = chol(sigma);
        if (flag)
            sigma = sigma + obj.PD_tol * eye(size(sigma)); 
        end

        REstimator = REstimator.InitConditions(mu, sigma);

        % UKF step
        [points, m, sigma] = REstimator.PropagationStep(step);

        % Plane propagation per particle
        plane = obj.PlanePropagation(last_epoch, new_epoch, [particles(1:pos-4,i); m(1:end-1,1)]);
        s = [plane(1:4,1); m];

        % Save the update
        particles(1:pos,i) = s;
        Sigma(pos-4:pos-1, pos-4:pos-1) = reshape(sigma, [pos/2 pos/2]);
        particles(pos+1:end,i) = reshape(Sigma, [], 1);
        aux = [repmat(plane(1:4,1), 1, size(points,2)); points];
        sigma_points(:,i) = reshape(aux, [], 1);
    end

    % Assemble the final propagated prior 
    PropPrior = [particles; weights];
end

%% Auxiliary functions 
% Dynamics 
function [state] = dynamics(epsilon, state, step)
    % Extract the (constant) plane state
    L = state(1,:);                   % Delaunay action
    G = state(2,:);                   % Angular momentum
    H = state(3,:);                   % Polar component of the angular momentum
    eta = G./L;                       % Eccentricity function

    % Propagate the anomalies
    omega = 1./L.^3 + (3/4) * epsilon./(L.^4.*eta.^3) .* (1-3*(H./G).^2);
    state(4,:) = state(4,:) + omega * step;
end