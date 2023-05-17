
function [PropPrior] = PropagationStep(obj, last_epoch, new_epoch, plane, Prior)
    % Constants 
    step = new_epoch - last_epoch;      % Time step
    step = step / obj.Tc;               % Dimensionalizing

    % Preallocation 
    particles = Prior(1:end-1,:);       % Prior states
    weights = Prior(end,:);             % Prior weights
    
    % Perform the propagation through the UKF step
    if (true && step > 0)
        % Extra the (constant) plane state
        L = plane(5);                   % Delaunay action
        G = plane(6);                   % Angular momentum
        H = plane(7);                   % Polar component of the angular momentum
        eta = G/L;                      % Eccentricity function

        % Propagate the anomalies
        omega = 1/L^3 + (3/4) * obj.epsilon/(L^4*eta^3) * (1-3*(H/G)^2);
        particles(1,:) = particles(1,:) + omega * step;

        % Propagate the covariance
        Sigma = 1/L^3 + (3/4) * obj.epsilon/(L^4*eta^3) * (1-3*(H/G)^2);
        particles(2,:) = particles(2,:) + Sigma * step;
    end

    % Assemble the final propagated prior 
    PropPrior = [repmat(plane,1,size(particles,2)); particles; weights];
end