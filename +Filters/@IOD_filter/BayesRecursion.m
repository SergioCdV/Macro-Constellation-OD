

function [f, X, N] = BayesRecursion(obj, tspan, Measurements)
    % Constants 

    % Preallocation 
    X = cell(1,length(tspan));
    N = cell(1,length(tspan)); 
    f = cell(1,length(tspan));

    % Bayesian recursion initialization
    last_epoch = tspan(1);                  % Initial epoch
    meas_index = 1;                         % Measurement index

    % State and weight initialization
    [particles, weights] = obj.Initialization();    
    tic

    % Main loop
    for i = 1:length(tspan)
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

        % Assemble the prior 
        Prior = [particles; weights];
        
        % Perform the correction step if new measurements are available 
        if (measurement_flag)
            for j = 1:new_measurements
                % Propagation step and weight proposal using the kinematic prior
                prop_epoch = Measurements{meas_index + j,1};
                if (last_epoch ~= prop_epoch)
                    [PropPrior] = obj.PropagationStep(last_epoch, prop_epoch, Prior);
                    PropPrior(end,:) = PropPrior(end,:) * obj.PS;
                else
                    PropPrior = Prior;
                end
                last_epoch = prop_epoch;

                % Transport the grid
                % [particles] = obj.TransportGrid(particles, X{i-1}(), obj.X);
    
                % Correction step 
                [Posterior] = obj.CorrectionStep(Measurements{meas_index + j,4}, Measurements{meas_index + j,5}, PropPrior);

                % Extraction
                particles = Posterior(1:end-1,:);
                weights = Posterior(end,:);
                weights = weights / sum(weights);
            end

            % Sanity check on the number of processed measurements 
            meas_index = meas_index + new_measurements;
        else
            % Propagate to the new epoch the clustered states
            prop_epoch = tspan(i);
            if (last_epoch ~= prop_epoch)
                % Particle propagation
                [Posterior] = obj.PropagationStep(last_epoch, prop_epoch, Prior);
                
                % Extraction
                particles = Posterior(1:end-1,:);
                weights = Posterior(end,:);
            else
                Posterior = Prior;
            end
            last_epoch = prop_epoch;

            % Transport the grid
            % [particles] = obj.TransportGrid(particles, X{i-1}(), obj.X);
        end

        % New pdf 
        f{i} = Posterior;
        toc

        % Estimation on the number of targets (number of planes in the constellation)
        obj.N = max(1, round( sum(weights) ));
        N{i} = obj.N;

        % State estimation
        obj.X = obj.StateEstimation(particles, weights, obj.N);
        X{i} = obj.X;

        % Posterior representation 
        f{i} = [particles; weights];

        % Check for resampling 
        Neff = 1/sum(weights.^2);

        if (Neff < round(2/3*size(particles,2)) || Neff > 1e5)
            M = (obj.L * obj.M + 1);
            particles = zeros(7, M * size(obj.X,2) ) ;
            for j = 1:size(obj.X,2)
                particles(1:4, 1+M*(j-1):M*j) = obj.UniformTangentQuat(obj.L, obj.M, obj.X(1:4,j));
                particles(5:7, 1+M*(j-1):M*j) = obj.AffineSampling(M, obj.X(5:7,j), reshape(obj.X(8:end,j), [3 3]));
            end

            % Physical constraint on the Delaunay action
            particles(5,:) = ones(1,size(particles,2)) .* (particles(5,:) <= 0) + particles(5,:) .* (particles(5,:) > 0);

            % Physical constraint on the eccentricity
            particles(6,:) = repmat(0.01, 1, size(particles,2)) .* (particles(2,:) < 0) + particles(2,:) .* (particles(1,:) >= 0);
        
            % Correlation
%             for j = 1:size(particles,2)
%                 STM = [1 0 0; ...
%                       sqrt(1-particles(6,j).^2) -particles(5,j)*particles(6,j)/(sqrt(1-particles(6,j).^2)) 0; ...
%                       sqrt(1-particles(6,j).^2) * particles(7,i) -particles(5,j)*particles(6,j)* particles(7,j)/(sqrt(1-particles(6,j).^2)) particles(5,j) .* sqrt(1-particles(6,j).^2)];
%                 
%                 particles(6,j) = particles(5,j) .* sqrt(1-particles(6,j).^2);
%                 particles(7,j) = particles(6,j) .* particles(7,j);
%                 particles(5:7,j) = particles(6:7,j) / det(STM);
%             end

            particles(6,:) = particles(5,:) .* sqrt(1-particles(6,:).^2);
            particles(7,:) = particles(6,:) .* particles(7,:);

            % New weights
            weights = repmat(1/size(particles,2), 1, size(particles,2));
        end
    end
end