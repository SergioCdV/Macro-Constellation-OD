
function [Posterior] = ExtendedCorrectionStep(obj, indices, Measurements, Estimator, PropPrior)
    % Create the measurements partitions and cells 
    if (length(indices) >= sqrt(obj.Gamma) * obj.Gamma * obj.N && length(indices) < 11)
        partitions_set = {};

        for i = length(indices):-1:3 * obj.N
            aux = SetPartition.SetPartition(indices,i);
            partitions_set = [partitions_set; aux];
        end

        partitions_set = [partitions_set; {num2cell(indices)}];

    elseif (length(indices) > 11)
        partitions_set = {num2cell(indices)};   
    elseif (length(indices) > 1)
        [partitions_set] = SetPartition.SetPartition(indices);
    else
        partitions_set = {num2cell(indices)};
    end

    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);

    pos = Estimator.StateDim+1;
    Estimator.StateDim = 4; 
    REstimator = Estimator.Init();
    
    % Corrector step of the KF for the covariance and means
    Sigma = reshape(particles(pos+1:pos+(pos-1)^2,1), [pos-1 pos-1]); 
    F = Sigma(1:3,1:3)^(-1);

    % Preallocation
    L = size(particles,2);                  % Number of particles
    DW = ones(1,size(partitions_set,1));
    update_particle = [];
    psi = [];
    psi_init = 1;
    psi_end = 0;

    % Summation over the partitions
    for i = 1:size(partitions_set,1)
        % Preallocation 
        p = partitions_set{i,:};        % Partition in the set

        % Arrange the macro-measurement and observation model
        for j = 1:size(p,2)
            W = p{j};                   % Cell in the partition
            W_card = size(W,2);         % Cardinality of the cell

            % Preallocation 
            Z = [];                                     % Measurement
            R = [];                                     % Covariance
            SensorModality = cell(W_card,1);            % Sensor modality
            ObservationModel = cell(W_card,1);          % Observation model
            Likelihood = cell(W_card,1);                % Likelihood function 
            l = zeros(1,W_card);                        % Likelihood evaluation

            for k = 1:W_card
                z = Measurements{W(k),2}(2:end).';                          % Individual measurement
                Z = [Z; z];                                                 % Measurement assemble

                r = reshape(Measurements{W(k),6}, [size(z,1) size(z,1)]);   % Individual covariance
                R = blkdiag(R, r);                                          % Covariance assemble

                % Functions 
                Likelihood{k} = Measurements{W(k),4};
                ObservationModel{k} = Measurements{W(k),5};
                SensorModality{k} = Measurements{W(k),7};
            end

            % Estimator 
            REstimator.R = R;
            REstimator = REstimator.AssignObservationProcess(size(Z,1), @(state)FullProcess(obj, SensorModality, ObservationModel, state));

            aux_particles = particles;
            aux_weights = weights;
            for k = 1:L
                Sigma = blkdiag(1, reshape( particles(pos+1:pos+(pos-1)^2,k), [pos-1 pos-1]));
                sigma_points = reshape( particles((pos-1)^2+pos+1:end,k), pos, []);
    
                % UKF step
                [mu, S, ~, ~] = REstimator.CorrectionStep(sigma_points, particles(1:pos,k), Sigma, Z);
                S = S(2:end,2:end);

                % ADMM step 
                X = projection_constraints(mu(5:7,1), box);
                
                % Particles update
                aux_particles(1:pos,index) = [mu(1:4,1); X; mu(end,:)];
                aux_particles(pos+1:pos+(pos-1)^2,k) = reshape(S, [], 1);

                % Likelihood update
                [y, dimensions] = FullProcess(obj, SensorModality, ObservationModel, mu);
                dim = 0;
                if (~isempty(y) && all(~isnan(y)))
                    for n = 1:W_card
                        % Evaluate the likelihood
                        l(n) = feval(Likelihood{n}, y(1+dim:dim+dimensions(n)));
                        dim = dim+dimensions(n);
                    end
                end

                % Weights update 
                Gamma = exp(-obj.Gamma) * obj.Gamma^W_card;
                aux_weights(k) = Gamma * obj.PD * prod(l) * weights(k);
            end

            dw = (1 * (W_card == 1) + sum(aux_weights));
            aux_weights = aux_weights / dw;

            DW(i) = DW(i) * dw;

            % Save the updated particles 
            update_particle = [update_particle aux_particles];
            psi = [psi aux_weights];
            psi_end = psi_end + L;
        end
        
        psi(1,psi_init:psi_end) = DW(i) * psi(1,psi_init:psi_end);
        psi_init = psi_end+1;
    end

    % Normalize 
    psi = psi / sum(DW,2);

    % Add the misdetected particles 
    psi = [weights psi];
    update_particle = [particles update_particle];

    % Update the weights with the non-detected particles 
    psi(1,1:L) = (1 - (1-exp(-obj.Gamma)) * obj.PD) .* weights(1,1:L);

    % Check for NaN 
    pos_nan = isnan(psi(1,:)); 
    psi(1,pos_nan) = zeros(1,sum(pos_nan));
    
    % Assemble the posterior 
    Posterior = [update_particle(1:pos+(pos-1)^2,:); psi];
end

%% Auxiliary functions 
% Full measurement process from the Delaunay set 
function [y, dimensions] = FullProcess(obj, SensorModality, ObservationModel, state)
    % Preallocation 
    for i = 1:size(state,2)

        measurement = []; 
        dimensions = [];

        for j = 1:size(SensorModality,1)
            X = projection_constraints(state(5:7,i), box);
            proj_state = [state(1:4,i); X; state(end,i)];
            State = obj.ParticleState(SensorModality{j}, proj_state).';
            [~, aux] = feval(ObservationModel{j}, State);

            % Dimensionalizing 
            switch (SensorModality{j})
                case 'RADAR'
                    aux = aux ./ [obj.Re obj.Re/obj.Tc];
                case 'INERTIAL'
                    aux = aux / obj.Re;
            end

            measurement = [measurement; aux.'];
            dimensions = [dimensions; length(aux)];
        end

        if (~isempty(measurement))
            y(1:length(measurement),i) = measurement;
        else
            y(1,i) = NaN;
        end
    end
end

% Projection step 
function [X] = projection_constraints(X, box)

    for k = 1:size(X,2)
        % Projection of the Delaunay action 
        if (X(1,k) < box(1,1))
            X(1,k) = box(1,1);
        elseif (X(1,k) > box(1,2))
            X(1,k) = box(1,2);
        end
        
        % Projection of the angular momentum 
        if (X(2,k) / X(1,k) < box(2,1))
            X(2,k) = X(1,k) * box(2,1);
        elseif (X(2,k) / X(1,k) > box(2,2))
            X(2,k) = X(1,k) * box(2,2);
        end
    
        % Projection of the nodal angular momentum 
        if (X(3,k) / X(2,k) < box(3,1))
            X(3,k) = X(2,k) * box(3,1);
        elseif (X(3,k) / X(2,k) > box(3,2))
            X(3,k) = X(2,k) * box(3,2);
        end
    end
end