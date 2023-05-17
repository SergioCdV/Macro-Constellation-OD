

function [f, X, N] = BayesRecursion(obj, tspan, Measurements)
    % Repeatibility 
    rng(1); 

    % Preallocation 
    X = cell(1,length(tspan));
    N = cell(1,length(tspan)); 
    f = cell(1,length(tspan));
    time = zeros(1,length(tspan));

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

            for j = 0:max(group)-1
                % Propagation step and weight proposal using the kinematic prior
                indices = index(group == j+1);
                prop_epoch = Measurements{meas_index + indices(1),1};
                if (last_epoch ~= prop_epoch)
                    [PropPrior] = obj.PropagationStep(last_epoch, prop_epoch, Prior);
                    PropPrior(end,:) = PropPrior(end,:) * obj.PS;
                else
                    PropPrior = Prior;
                end
                last_epoch = prop_epoch;

                % Birth particles
                born_particles = obj.Birth();
                PropPrior(:,size(particles,2)+1:size(particles,2)+size(born_particles,2)) = born_particles;
    
                % Correction step 
                [Posterior] = obj.CorrectionStep(Measurements, PropPrior, meas_index+indices);
            end

            % Sanity check on the number of processed measurements 
            meas_index = meas_index + new_measurements;
        else
            % Propagate to the new epoch the clustered states
            prop_epoch = tspan(i);
            if (last_epoch ~= prop_epoch)
                % Particle propagation
                [Posterior] = obj.PropagationStep(last_epoch, prop_epoch, Prior);
            else
                Posterior = Prior;
            end
            last_epoch = prop_epoch;
        end

        % New pdf 
        f{i} = Posterior;
        particles = Posterior(1:end-1,:);
        weights = Posterior(end,:);

        % Estimation on the number of targets (number of spacecraft in the planes)
        T = sum(weights);
        obj.N = round(T);
        
        obj.X = obj.StateEstimation(particles, weights, T);
        X{i} = obj.X;
        N{i} = obj.N;

        % Kalman update for the orbital plane states

        % New prior 
        Prior = [particles; weights];

        time(i) = toc;
        fprintf('Iteration running time: %.4f s\n', time(i));

    end
    fprintf('------------------------------\n');
    fprintf('Bayes filter recursion finished.\n');
    fprintf('Total running time: %.4f s\n', sum(time));
end