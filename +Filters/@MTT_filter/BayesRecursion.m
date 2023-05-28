

function [f, X, N, Prior] = BayesRecursion(obj, tspan, Measurements)
    % Repeatibility 
    rng(1); 

    % Preallocation 
    X = cell(1,length(tspan));
    N = cell(1,length(tspan)); 
    f = cell(1,length(tspan));
    time = zeros(1,length(tspan));

    % Preallocate the estimator 
    pos = 8;
    Q = 1e-7 * eye(pos);

    switch (obj.KF_type)
        case 'EKF'
            AnomalyEstimator = Filters.EKF();
        case 'UKF-A'
            AnomalyEstimator = Filters.UKF('UKF-A', 2, 1E-1, 0);
            AnomalyEstimator.StateDim = pos-1;
            AnomalyEstimator = AnomalyEstimator.AdditiveCovariances(Q, zeros(3)).Init();
        case 'UKF-S'
            AnomalyEstimator = Filters.UKF('UKF-S', 2, 1E-1, 0);
            AnomalyEstimator.StateDim = pos-1;
            AnomalyEstimator = AnomalyEstimator.AdditiveCovariances(Q, zeros(3)).Init();
    end

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

                % Augmentation of the perifocal attitude particle representation via tangent space mapping
                [Prior, obj.Gibbs_vector] = obj.PerifocalQuatSampling(Prior);

                [Prior(1:pos-1,:)] = obj.PlanePropagation(last_epoch, prop_epoch, Prior(1:pos-1,:));
                [PropPrior, sigma_points] = obj.PropagationStep(last_epoch, prop_epoch, AnomalyEstimator, Prior);
                PropPrior = [PropPrior(1:end-1,:); sigma_points; PropPrior(end,:)];
                PropPrior(end,:) =  obj.PS * PropPrior(end,:);

                % Birth particles
                born_particles = obj.Birth();
                [born_particles, bsigma_points] = obj.PropagationStep(last_epoch, last_epoch, AnomalyEstimator, born_particles);
                born_particles = [born_particles(1:end-1,:); bsigma_points; born_particles(end,:)];

                PropPrior = [PropPrior born_particles];
    
                % Correction step 
                [Posterior] = obj.CorrectionStep(meas_index+indices, Measurements, AnomalyEstimator, PropPrior);
%                 [Posterior] = obj.ExtendedCorrectionStep(meas_index+indices, Measurements, AnomalyEstimator, PropPrior);
                Prior = Posterior;
                last_epoch = prop_epoch;
            end

            % Final posterior
            particles = Posterior(1:end-1,:);
            weights = Posterior(end,:);

            % Averaging of the perifocal frame (ML state estimation)
            [plane, Sigma] = obj.PerifocalUpdate(weights / sum(weights), particles(1:pos-1,:));

            obj.planes = [plane; reshape(Sigma, [], 1)];
            particles(1:pos-1,:) = repmat(obj.planes(1:pos-1,:), 1, size(particles,2));
            
            for j = 1:size(particles,2)
                Sigma_t = reshape(particles(pos+1:end,j), [pos-1 pos-1]);
                Sigma_t(1:end-1,1:end-1) = Sigma;
                particles(pos+1:end,j) = reshape(Sigma_t, [], 1);
            end
            
            % Estimation on the number of targets per plane
            T = sum(weights,2);
            obj.N = round( T + (1-obj.PD) * obj.PS * obj.N );

            % Sanity check on the number of processed measurements 
            meas_index = meas_index + new_measurements;

            % Pruning of the weights and merging of the particles
            [particles, weights] = obj.Pruning(particles, weights); 

        else
            % Propagate to the new epoch the clustered states
            prop_epoch = tspan(i);
            if (last_epoch ~= prop_epoch)
                % Plane propagation 
                [Prior(1:pos-1,1)] = obj.PlanePropagation(last_epoch, prop_epoch, Prior(1:pos-1,1));
                Prior(1:pos-1,2:end) = repmat(Prior(1:pos-1,1), 1, size(Prior,2)-1);

                % Particle propagation
                [Posterior] = obj.PropagationStep(last_epoch, prop_epoch, AnomalyEstimator, Prior);
            else
                Posterior = Prior;
            end

            particles = Posterior(1:end-1,:);
            weights = Posterior(end,:);

            last_epoch = prop_epoch;
        end

        future = i + max(1,new_measurements);

        % New pdf 
        aux = zeros(length(obj.nu),size(particles,2));
        for j = 1:size(particles,2)
            aux(:,j) = obj.wrapped_normal(1e-7, obj.nu.', mod(particles(pos,j),2*pi), particles(end,j));
        end
        f{i} = sum(weights .* aux,2);

        % Estimation on the anomaly space
        obj.X = obj.StateEstimation(particles, weights, obj.N);
        X{i} = obj.X;
        N{i} = obj.N;

        for j = 1:new_measurements-1
            f{i+j} = f{i};
            X{i+j} = obj.X;
            N{i+j} = obj.N;
        end

        % New prior 
        Prior = [particles; weights];

        time(i) = toc;
        fprintf('Iteration running time: %.4f s\n', time(i));

        % Update the timer 
        i = future;
    end

    fprintf('------------------------------\n');
    fprintf('Bayes filter recursion finished.\n');
    fprintf('Total running time: %.4f s\n', sum(time));
end