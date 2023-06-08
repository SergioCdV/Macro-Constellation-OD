
function [PropPrior, sigma_points] = PropagationStep(obj, last_epoch, new_epoch, PlaneEstimator, Prior)
    % Constants 
    step = new_epoch - last_epoch;      % Time step
    step = 86400 / obj.Tc * step;       % Dimensionalizing
    pos = 7;                            % State dimension
    particles = Prior(1:end-1,:);       % Prior states
    weights = Prior(end,:);             % Prior weights

    % Preallocation 
    aug_particles = zeros(pos+3+(pos-1)^2, size(particles,2));
    sigma_points = zeros((pos+3) * (2 * PlaneEstimator.StateDim + 1), size(particles,2));

    % Perform the propagation through the KF step
    for i = 1:size(particles,2)
        % Initial conditions
        mu = particles(1:pos,i);
        Sigma = reshape(particles(pos+1:end,i), [pos-1 pos-1]); 

        PlaneEstimator = PlaneEstimator.InitConditions(mu, Sigma);

        % UKF step
        [points, m, sigma] = PlaneEstimator.PropagationStep(step);

        % Save the update
        Sigma = reshape(sigma, [pos-1 pos-1]);
        aug_particles(1:pos+3,i) = m;
        aug_particles(pos+3+1:end,i) = reshape(Sigma, [], 1);
        sigma_points(:,i) = reshape(points, [], 1);
    end

    % Assemble the final propagated prior 
    PropPrior = [aug_particles; weights];
end
