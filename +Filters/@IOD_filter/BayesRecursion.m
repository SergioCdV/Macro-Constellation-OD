

function [f, X, N] = BayesRecursion(obj, tspan, Measurements)
    % Repeatibility 
    rng(1); 

    % Preallocation 
    X = cell(1,length(tspan));
    N = cell(1,length(tspan)); 
    f = cell(1,length(tspan));

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
    for i = 1:100%length(tspan)
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

                % Transport the grid
                % [particles] = obj.TransportGrid(particles, X{i-1}(), obj.X);

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

            % Transport the grid
            % [particles] = obj.TransportGrid(particles, X{i-1}(), obj.X);
        end

        % New pdf 
        f{i} = Posterior;
        particles = Posterior(1:end-1,:);
        weights = Posterior(end,:);

        % Estimation on the number of targets (number of planes in the constellation)
        T = sum(weights);
        obj.N = round(T);
        N{i} = obj.N;

        % Check for resampling 
        Neff = 1/sum(weights.^2);

        if (Neff < (2/3) * size(particles,2 || size(particles,2) > obj.Jmax))        
            % Resampling 
%             if (size(particles,2) > obj.Jmax)
%                 % Pruning 
%                 % [particles, weights] = Pruning(obj, particles, weights);
% 
%                 % State estimation
% %             if (obj.N)
% %                 obj.X = obj.StateEstimation(particles, weights, obj.N);
% %                 X{i} = obj.X;
% %             end
% 
%             else
                % State estimation
                if (obj.N)
                    obj.X = obj.StateEstimation(particles, weights, obj.N);
                    X{i} = obj.X;
                end

                % Resampling
                J = max(obj.N,1) * (1 + obj.L * obj.M);
                [particles, weights] = obj.Resampling(particles, weights/T, J);
                weights = weights * T;    
%             end   

        elseif (obj.N)
            % State estimation
            obj.X = obj.StateEstimation(particles, weights, obj.N);
            X{i} = obj.X;
    
            % Double-covering of the sphere
            % obj.X = [obj.X [-obj.X(1:4,:); obj.X(5:end,:)]];

            M = (obj.L * obj.M + 1);
            particles = zeros(7, M * size(obj.X,2) ) ;
            for j = 1:size(obj.X,2)
                particles(1:4, 1+M*(j-1):M*j) = obj.UniformTangentQuat(obj.L, obj.M, obj.X(1:4,j));
                actions = obj.GibbsSampling(2 * M, obj.X(5:7,j), reshape(obj.X(8:end,j), [3 3]), [1 7; 0.9 6.9; -7 7]);
                particles(5:7, 1+M*(j-1):M*j) = actions(:,M+1:end);
            end

            % New weights
            weights = repmat(T/size(particles,2), 1, size(particles,2));
        end

        % New prior 
        Prior = [particles; weights];

        fprintf('Iteration running time: %.4f s\n', toc);

    end
    fprintf('------------------------------\n');
    fprintf('Recursion finished.\n');
end