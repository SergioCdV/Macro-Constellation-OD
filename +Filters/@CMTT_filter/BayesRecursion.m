

function [f, X, N, Prior, E] = BayesRecursion(obj, tspan, Measurements)
    % Repeatibility 
    rng(1); 

    % Preallocation 
    X = cell(1,length(tspan));
    N = zeros(1,length(tspan)); 
    f = cell(1,length(tspan));
    time = zeros(1,length(tspan));
    E = zeros(1,length(tspan));

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

    Q = 1e-7 * eye(pos-2);
    switch (obj.PKF_type)
        case 'UKF-A'
            PlaneEstimator = Filters.USQUE('UKF-A', 2, 1E-1, 0, 1);
            PlaneEstimator.StateDim = pos-2;
            PlaneEstimator = PlaneEstimator.AdditiveCovariances(Q, zeros(3)).Init();
        case 'UKF-S'
            PlaneEstimator = Filters.USQUE('UKF-S', 2, 1E-1, 0, 1);
            PlaneEstimator.StateDim = pos-2;
            PlaneEstimator = PlaneEstimator.AdditiveCovariances(Q, zeros(3)).Init();
    end

    % Bayesian recursion initialization
    last_epoch = tspan(1);                  % Initial epoch
    meas_index = 1;                         % Measurement index

    fprintf('Initialization... \n');
    fprintf('------------------------------\n');
    fprintf('Running the filter: \n');

    % State and weight initialization
    [AnomalyParticles, AnomalyWeights] = obj.Initialization();  

    % Assemble the prior 
    Prior = [AnomalyParticles; AnomalyWeights];

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

                % Plane propagation 
                [PropPlane, PropSigmaPlane, plane_points] = obj.PlanePropagation(last_epoch, prop_epoch, PlaneEstimator, obj.planes);
                Prior(1:pos-1,:) = repmat(PropPlane([pos-1:pos+2 4:6],1), 1, size(Prior,2));
                obj.Gibbs_vector = repmat(PropPlane(1:3,1), 1, size(Prior,2));
                for j = 1:size(AnomalyParticles,2)
                    Sigma_t = reshape(AnomalyParticles(pos+1:end,j), [pos-1 pos-1]);
                    Sigma_t(1:end-1,1:end-1) = PropSigmaPlane;
                    AnomalyParticles(pos+1:end,j) = reshape(Sigma_t, [], 1);
                end

                % Anomaly propagation
                [PropPrior, sigma_points] = obj.PropagationStep(last_epoch, prop_epoch, AnomalyEstimator, Prior);
                PropPrior = [PropPrior(1:end-1,:); sigma_points; PropPrior(end,:)];
                PropPrior(end,:) = obj.PS * PropPrior(end,:);

                % Birth particles
                born_particles = obj.Birth();
                [born_particles, bsigma_points] = obj.PropagationStep(last_epoch, last_epoch, AnomalyEstimator, born_particles);
                born_particles = [born_particles(1:end-1,:); bsigma_points; born_particles(end,:)];

                PropPrior = [PropPrior born_particles];
    
                % Correction step 
                if (obj.ExtendedTarget)
                    [Posterior] = obj.ExtendedCorrectionStep(meas_index+indices, Measurements, AnomalyEstimator, PropPrior);
                else
                    [Posterior] = obj.CorrectionStep(meas_index+indices, Measurements, AnomalyEstimator, PropPrior);
                end

                % Estimation of the perifocal frame
                AnomalyParticles = Posterior(1:end-1,:);
                AnomalyWeights = Posterior(end,:);
                [plane, Sigma] = obj.PlaneCorrection(PlaneEstimator, PropPlane, PropSigmaPlane, plane_points, AnomalyParticles);
    
                obj.planes = [plane; reshape(Sigma, [], 1)];
                AnomalyParticles(1:pos-1,:) = repmat(obj.planes(1:pos-1,:), 1, size(AnomalyParticles,2));
                
                for j = 1:size(AnomalyParticles,2)
                    Sigma_t = reshape(AnomalyParticles(pos+1:end,j), [pos-1 pos-1]);
                    Sigma_t(1:end-1,1:end-1) = Sigma;
                    AnomalyParticles(pos+1:end,j) = reshape(Sigma_t, [], 1);
                end

                Prior = Posterior;
                last_epoch = prop_epoch;
            end
            
            % Estimation on the number of targets per plane
            T = sum(AnomalyWeights,2);
            obj.N = round( T + (1-obj.PD) * obj.PS * obj.N );

            % Sanity check on the number of processed measurements 
            meas_index = meas_index + new_measurements;

            % Pruning of the weights and merging of the particles
            [AnomalyParticles, AnomalyWeights] = obj.Pruning(AnomalyParticles, AnomalyWeights); 

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

            AnomalyParticles = Posterior(1:end-1,:);
            AnomalyWeights = Posterior(end,:);

            last_epoch = prop_epoch;
        end

        future = i + max(1,new_measurements);

        % New pdf 
        aux = zeros(length(obj.nu),size(AnomalyParticles,2));
        for j = 1:size(AnomalyParticles,2)
            aux(:,j) = obj.wrapped_normal(1e-7, obj.nu.', mod(AnomalyParticles(pos,j),2*pi), AnomalyParticles(end,j));

            % Entropy characterization 
            Sigma = reshape(AnomalyParticles(pos+1:end,j), [pos-1 pos-1]);
            entropy = 0.5 * log(det(2*pi*exp(1)*Sigma));
            E(i) = entropy * (entropy < E(i)) + E(i) * (entropy >= E(i));
        end
        f{i} = sum(AnomalyWeights .* aux,2);

        % Estimation on the anomaly space
        obj.X = obj.StateEstimation(AnomalyParticles, AnomalyWeights, obj.N);
        X{i} = obj.X;
        N(i) = obj.N;

        for j = 1:new_measurements-1
            f{i+j} = f{i};
            X{i+j} = obj.X;
            N(i+j) = obj.N;
            E(i+j) = E(i);
        end

        % New prior 
        Prior = [AnomalyParticles; AnomalyWeights];

        time(i) = toc;
        fprintf('Iteration running time: %.4f s\n', time(i));

        % Update the timer 
        i = future;
    end

    fprintf('------------------------------\n');
    fprintf('Bayes filter recursion finished.\n');
    fprintf('Total running time: %.4f s\n', sum(time));
end