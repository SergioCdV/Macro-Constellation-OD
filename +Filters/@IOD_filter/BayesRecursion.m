

function [f, X, N] = BayesRecursion(obj, tspan, Measurements)
    % Repeatibility 
    rng(1); 

    % Preallocation 
    X = cell(1,length(tspan));
    N = cell(1,length(tspan)); 
    f = cell(1,length(tspan));
    time = zeros(1,length(tspan));

    % Quadrature precomputation 
%     CC_quad = Numerics.CollocationMesh.ChebyshevGrid(1e2);
%     [obj.dt(1,:), obj.dt(2,:)] = CC_quad.Domain(0, 2*pi, CC_quad.tau);
%     obj.W = CC_quad.W;
%     obj.nu = CC_quad.tau;

    obj.nu = linspace(0, 2*pi, 3e2);

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

            for j = 0:max(group)-1
                % Propagation step and weight proposal using the kinematic prior
                indices = index(group == j+1);
                prop_epoch = Measurements{meas_index + indices(1),1};

                [PropPrior] = obj.PropagationStep(last_epoch, prop_epoch, Prior);
                PropPrior(end,:) =  obj.PS * PropPrior(end,:);

                % Birth particles
                born_particles = obj.Birth();
                PropPrior = [born_particles PropPrior];
    
                % Correction step 
                [Posterior] = obj.CorrectionStep(meas_index+indices, Measurements, PropPrior);

                Prior = Posterior;
                last_epoch = prop_epoch;
            end

            % Final posterior
            particles = Posterior(1:end-1,:);
            weights = Posterior(end,:);

            % Estimation on the number of targets (number of planes in the constellation)
            T = sum(weights,2);

            % Estimation of the plane perifocal state
            obj.X = obj.StateEstimation(particles, weights, T);
            obj.N = size(obj.X,2);

            % Sanity check on the number of processed measurements 
            meas_index = meas_index + new_measurements;
            
        else
            % Propagate to the new epoch the clustered states
            if (last_epoch ~= prop_epoch)
                % Particle propagation
                [Posterior] = obj.PropagationStep(last_epoch, prop_epoch, Prior);
            else
                Posterior = Prior;
            end

            particles = Posterior(1:end-1,:);
            weights = Posterior(end,:);

            last_epoch = prop_epoch;
        end

        future = i + max(1,new_measurements);

        % New pdf 
        f{i} = Posterior;

        % Resampling the actions
        M = (1 + obj.L * obj.M);
        J = max(obj.N,1) * M;
        [particles, weights] = obj.Resampling(particles, weights/T, J);

        Neff = 1/sum(weights.^2,2);

        if (Neff < (1/3) * size(particles,2))
            % Re-population of the quaternions
            for j = 1:size(obj.X,2)
                [qd, a] = obj.PerifocalQuatSampling(obj.X(:,j));
                particles(:, 1+M*(j-1):M*j) = repmat(qd(1:7,:), 1, M);

                % Conditioning
                Sigma = reshape(obj.X(8:end,j), [6 6]);
                F = Sigma(1:3,1:3)^(-1);
                sigma = Sigma(4:6,4:6) - Sigma(4:6,1:3) * F * Sigma(4:6,1:3).';
                mu = obj.X(5:7,j) + Sigma(4:6,1:3) * F * a;

                % Gaussian sampling to improve impovershment in the action space
                actions = obj.GibbsSampling(2 * M, mu, sigma, obj.search_limit);
                % actions = mvnrnd(mu, sigma, M).'; 
                particles(5:7, 1+M*(j-1):M*j) = actions(:,M+1:end);
            end
        end

        % Normalizing the weights 
        weights = weights * T;

        % Save results
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