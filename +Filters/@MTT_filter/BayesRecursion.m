

function [f, X, N] = BayesRecursion(obj, tspan, Measurements)
    % Repeatibility 
    rng(1); 

    % Preallocation 
    X = cell(1,length(tspan));
    N = cell(1,length(tspan)); 
    f = cell(1,length(tspan));
    time = zeros(1,length(tspan));

    % Preallocate the estimator 
    Q = zeros(8);
    Q(end,end) = 1e-3;
    switch (obj.KF_type)
        case 'EKF'
            AnomalyEstimator = Filters.EKF();
        case 'UKF-A'
            AnomalyEstimator = Filters.UKF('UKF-A', 2, 1E-2, 0);
            AnomalyEstimator.StateDim = 8;
            AnomalyEstimator = AnomalyEstimator.AdditiveCovariances(Q, zeros(3)).Init();
        case 'UKF-S'
            AnomalyEstimator = Filters.UKF('UKF-S', 2, 1E-2, 0);
            AnomalyEstimator.StateDim = 8;
            AnomalyEstimator = AnomalyEstimator.AdditiveCovariances(Q, zeros(3)).Init();
    end

    % Quadrature precomputation 
    obj.nu = linspace(0,2*pi,1e2);          % Anomaly space

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
    for i = 1:length(tspan)
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
            % Group the measurement per their epoch
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
                [PropPrior, sigma_points] = obj.PropagationStep(last_epoch, prop_epoch, AnomalyEstimator, Prior);
                PropPrior = [PropPrior(1:end-1,:); sigma_points; PropPrior(end,:)];
                PropPrior(end,:) = PropPrior(end,:) * obj.PS;
                last_epoch = prop_epoch;

                % Birth particles
                born_particles = obj.Birth();
                [born_particles, bsigma_points] = obj.PropagationStep(last_epoch, last_epoch, AnomalyEstimator, born_particles);
                born_particles = [born_particles(1:end-1,:); bsigma_points; born_particles(end,:)];
                PropPrior = [PropPrior born_particles];
    
                % Correction step 
                [Posterior] = obj.CorrectionStep(meas_index+indices, Measurements, AnomalyEstimator, PropPrior);
            end

            % Sanity check on the number of processed measurements 
            meas_index = meas_index + new_measurements;
        else
            % Propagate to the new epoch the clustered states
            prop_epoch = tspan(i);
            if (last_epoch ~= prop_epoch)
                % Particle propagation
                [Posterior] = obj.PropagationStep(last_epoch, prop_epoch, AnomalyEstimator, Prior);
            else
                Posterior = Prior;
            end
            last_epoch = prop_epoch;
        end

        % Kalman update for the orbital plane states
        if (measurement_flag)
        
        end

        % New pdf 
        f{i} = Posterior;
        particles = Posterior(1:end-1,:);
        weights = Posterior(end,:);

        % Estimation on the number of targets per plane
        T = sum(weights);
        obj.N = round(T);
        
        [particles, weights] = obj.Pruning(particles, weights);            % Prunning

        obj.X = obj.StateEstimation(particles, weights, obj.N);
        X{i} = obj.X;
        N{i} = obj.N;

        % New prior 
        Prior = [particles; weights];

        time(i) = toc;
        fprintf('Iteration running time: %.4f s\n', time(i));

    end
    fprintf('------------------------------\n');
    fprintf('Bayes filter recursion finished.\n');
    fprintf('Total running time: %.4f s\n', sum(time));
end