

function [X, N, Prior, E] = BayesRecursion(obj, tspan, Measurements)
    % Repeatibility 
    rng(1); 

    % Preallocation 
    X = cell(1,length(tspan));          % State estimation
    N = zeros(1,length(tspan));         % Number of planes
    time = zeros(1,length(tspan));      % Time for each iteration
    E = zeros(1,length(tspan));         % Differential entropy

    % Preallocate the estimator 
    pos = 7;
    Q = 1e-4 * eye(pos-1);

    switch (obj.KF_type)
        case 'UKF-A'
            PlaneEstimator = Filters.USQUE('UKF-A', 2, 1E-1, 0, 1);
            PlaneEstimator.StateDim = pos-1;
            PlaneEstimator = PlaneEstimator.AdditiveCovariances(Q, zeros(3)).Init();
        case 'UKF-S'
            PlaneEstimator = Filters.USQUE('UKF-S', 2, 1E-1, 0, 1);
            PlaneEstimator.StateDim = pos-1;
            PlaneEstimator = PlaneEstimator.AdditiveCovariances(Q, zeros(3)).Init();
    end

    PlaneEstimator = PlaneEstimator.AssignStateProcess(PlaneEstimator.StateDim, @(planes, step)PlanePropagation(obj, planes, step)).Init();

    % Bayesian recursion initialization
    last_epoch = tspan(1);                  % Initial epoch
    meas_index = 1;                         % Measurement index

    fprintf('Initialization... \n');
    fprintf('------------------------------\n');
    fprintf('Running the filter: \n');

    % State and weight initialization
    [particles, weights] = obj.Initialization();  

    % Assemble the prior 
    Prior = [particles; weights];

    % Main loop
    T = [0; obj.N];
    i  = 1;
    while (i <= length(tspan))
        tic 

        % Check for new measurements to process
        new_measurements = 0;

        if ( ~isempty(Measurements) )
            GoOn = true;
            while (GoOn)
                if (meas_index + new_measurements > size(Measurements,1))
                    new_measurements = size(Measurements,1) - meas_index;
                    GoOn = false;
                elseif (Measurements{meas_index + new_measurements,1} <= tspan(i))
                    new_measurements = new_measurements + 1; 
                else
                    GoOn = false;
                end
            end
    
            % Enable the corrector step
            if ( new_measurements > 0 ) 
                measurement_flag = true;
            else
                measurement_flag = false;
            end
        else
            measurement_flag = false;
        end

        % Perform the correction step if new measurements are available 
        if (measurement_flag)
            % Group the measurement by their epoch
            group = zeros(1,new_measurements);

            index = 0:new_measurements-1;
            j = 1;
            while (~isempty(index) && (j <= new_measurements))
                cur_pos = index(1);
                res_indices = index(2:end);
                for k = 1:length(res_indices)
                    if (Measurements{meas_index + cur_pos,1} == Measurements{meas_index + res_indices(k),1})
                        group(res_indices(k) + 1) = j;
                        group(cur_pos + 1) = j;
                        index = index(group(index + 1) ~= j);
                    end
                end

                if (group(cur_pos + 1) == 0)
                    group(cur_pos + 1) = j;
                    index = index(group(index + 1) ~= j);
                end
                j = j+1;
            end

            index = 0:new_measurements-1;

            % Anomaly distribution update
            for j = 0:max(group)-1
                % Propagation step and weight proposal using the kinematic prior
                indices = index(group == j+1);
                prop_epoch = Measurements{meas_index + indices(1),1};

                [PropPrior, sigma_points] = obj.PropagationStep(last_epoch, prop_epoch, PlaneEstimator, Prior);
                PropPrior = [PropPrior(1:end-1,:); sigma_points; PropPrior(end,:)];
                PropPrior(end,:) =  obj.PS * PropPrior(end,:);

                % Birth particles
                born_particles = obj.Birth();
                [born_particles, bsigma_points] = obj.PropagationStep(last_epoch, last_epoch, PlaneEstimator, born_particles);
                born_particles = [born_particles(1:end-1,:); bsigma_points; born_particles(end,:)];

                PropPrior = [PropPrior born_particles];
    
                % Correction step 
                [Posterior] = obj.CorrectionStep(meas_index+indices, Measurements, PlaneEstimator, PropPrior);

                Prior = Posterior;
                last_epoch = prop_epoch;
            end

            % Final posterior
            particles = Posterior(1:end-1,:);
            weights = Posterior(end,:);
            
            % Estimation on the number of targets per plane
            T(1) = sum(weights,2);
            T(2) = round( T(1) + (1-obj.PD) * obj.PS * T(2) );

            % Sanity check on the number of processed measurements 
            meas_index = meas_index + new_measurements;

            % Estimation on the plane space
            obj.X = obj.StateEstimation(particles, weights, T(1));

        else
            % Propagate to the new epoch the clustered states
            prop_epoch = tspan(i);
            if (last_epoch ~= prop_epoch)
                % Particle propagation
                [Posterior] = obj.PropagationStep(last_epoch, prop_epoch, PlaneEstimator, Prior);
                Posterior = [Posterior(pos:pos+3,:); Posterior(4:pos-1,:); Posterior(pos+3+1:pos+3+(pos-1)^2,:); Posterior(end,:)];
            else
                Posterior = Prior;
            end
            
            particles = Posterior(1:end-1,:);
            weights = Posterior(end,:);

            last_epoch = prop_epoch;
        end

        future = i + max(1,new_measurements);

        % New pdf 
        for j = 1:size(particles,2)
            % Entropy characterization 
            Sigma = reshape(particles(pos+1:end,j), [pos-1 pos-1]);
            entropy = 0.5 * log(det(2*pi*exp(1)*Sigma));
            E(i) = entropy * (entropy < E(i)) + E(i) * (entropy >= E(i));
        end

        X{i} = obj.X;
        obj.N = size(obj.X,2);
        N(i) = obj.N;
        for j = 1:new_measurements-1
            X{i+j} = obj.X;
            N(i+j) = obj.N;
            E(i+j) = E(i);
        end

        % New prior 
        M = min(obj.M * obj.N, obj.Jmax);
        [particles, weights] = obj.Resampling(particles, weights / T(1), M);          % Pruning of the weights and merging of the particles via resampling
        Prior = [particles; T(1) * weights];

        time(i) = toc;
        fprintf('Iteration running time: %.4f s\n', time(i));

        % Update the timer 
        i = future;
    end

    fprintf('------------------------------\n');
    fprintf('Bayes filter recursion finished.\n');
    fprintf('Total running time: %.4f s\n', sum(time));
end