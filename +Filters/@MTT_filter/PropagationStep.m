
function [PropPrior, sigma_points] = PropagationStep(obj, last_epoch, new_epoch, Estimator, Prior)
    % Constants 
    step = new_epoch - last_epoch;      % Time step
    step = step / obj.Tc;               % Dimensionalizing
    pos = 8; 

    % Preallocation 
    particles = Prior(1:end-1,:);       % Prior states
    weights = Prior(end,:);             % Prior weights
    
    % Reduce the dimensionality dynamics 
    REstimator = Estimator;
    REstimator.StateDim = 4; 
    REstimator.Q = REstimator.Q([pos-3:pos], [pos-3:pos]);

    % Perform the propagation through the KF step
    REstimator = REstimator.Init().AssignStateProcess(pos, @(M, step)dynamics(obj.epsilon, M, step));
    sigma_points = zeros(size(particles));
    for i = 1:size(particles,2)
        % KF step 
        Sigma = reshape(particles(pos+1:end,i), [pos pos]); 
        sigma = Sigma([pos-3:pos], [pos-3:pos]);
        REstimator = REstimator.InitConditions(particles(pos-3:pos,i), sigma);
        [points, m, sigma] = REstimator.PropagationStep(step);
        Sigma([pos-3:pos], [pos-3:pos]) = reshape(sigma, [pos/2 pos/2]);
        particles(pos-3:pos,i) = m;
        particles(pos+1:pos+pos^2,i) = reshape(Sigma, [], 1);
        aux = [repmat(particles(1:pos-4,i), 1, size(points,2)); points];
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