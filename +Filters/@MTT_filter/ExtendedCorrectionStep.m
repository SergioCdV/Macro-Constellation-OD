
function [Posterior] = ExtendedCorrectionStep(obj, indices, Measurements, Estimator, PropPrior)
    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);

    pos = Estimator.StateDim+1;
    Estimator.StateDim = 4; 
    REstimator = Estimator.Init();
    
    % Corrector step of the KF for the covariance and means
    Sigma = reshape(particles(pos+1:pos+(pos-1)^2,1), [pos-1 pos-1]); 
    F = Sigma(1:3,1:3)^(-1);

    % Create the measurements partitions and cells 
    [partitions_set] = SetPartition.SetPartition(indices);

    % Preallocation
    L = size(particles,2);
    update_particle = [];
    psi = [];
    DW = ones(1,size(partitions_set,1));

    % Summation over the partitions
    for i = 1:size(partitions_set,1)
        % Preallocation 
        p = partitions_set{i,:};        % Partition in the set

        % Arrange the macro-measurement and observation model
        for j = 1:size(p,2)
            W = p{j};                   % Cell in the partition
            W_card = size(W,2);         % Cardinality of the cell

            % Preallocation 
            Z = [];
            R = [];
            SensorModality = cell(W_card,1);
            ObservationModel = cell(W_card,1);
            Likelihood = cell(W_card,1);

            for k = 1:W_card
                z = Measurements{W(k),2}(2:end).';                          % Individual measurement
                Z = [Z; z];                                                 % Measurement assemble

                r = reshape(Measurements{W(k),6}, [size(z,1) size(z,1)]);   % Individual covariance
                R = blkdiag(R, r);                                          % Covariance assemble

                % Functions 
                Likelihood{i} = Measurements{W(k),4};
                ObservationModel{i} = Measurements{W(k),5};
                SensorModality{i} = Measurements{W(k),7};
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
                [mu, S, ~, y] = REstimator.CorrectionStep(sigma_points, particles(1:pos,k), Sigma, Z);
                Sigma = Sigma(2:end,2:end);
                S = S(2:end,2:end);
    
                % Deconditioning
                mu(5:end) = mu(5:end) - Sigma(4:end,1:3) * F * obj.Gibbs_vector(:,i);
                S(4:end,4:end) = S(4:end,4:end) + S(4:end,1:3) * F * S(4:end,1:3).';

                % Particles update
                aux_particles(1:pos,k) = mu;
                aux_particles(pos+1:pos+(pos-1)^2,k) = reshape(S, [], 1);

                % Likelihood update
                if (isempty(y) || any(isnan(y)))
                    l = 0; 
                else
                    l = feval(Likelihood, y) * l;
                end

                % Weights update 
                Gamma = exp(-obj.Gamma) * obj.Gamma^W_card;
                aux_weights(k) = Gamma * obj.PD * l * aux_weights(k);
            end

            dw = (1 * (W_card == 1) + sum(aux_weights));
            aux_weights = aux_weights / dw;

            DW(i) = DW(i) * dw;

            % Save the updated particles 
            update_particle = [update_particle aux_particles];
            psi = [psi aux_weights];
        end
        
        psi(1,1+L*size(p,2)*(i-1):L*size(p,2)*i) = DW(i) * psi(1,1+L*size(p,2)*(i-1):L*size(p,2)*i);
    end

    % Check for NaN 
    psi = psi / sum(DW,2); 
    psi = [weights psi];
    update_particle = [particles update_particle];
    pos_nan = isnan(psi(1,:)); 
    psi(1,pos_nan) = zeros(1,sum(pos_nan));

    % Update the weights with the non-detected particles 
    psi(1,1:L) = (1 - (1-exp(-obj.gamma)) * obj.PD) .* weights(1,:);
    
    % Assemble the posterior 
    Posterior = [update_particle; psi];
end

%% Auxiliary functions 
% Full measurement process from the Delaunay set 
function [y] = FullProcess(obj, SensorModality, ObservationModel, state)
    % Preallocation 
    for i = 1:size(state,2)
        State = obj.ParticleState(state(:,i)).';
        [~, aux] = feval(ObservationModel, State);
        if (~isempty(aux))
            % Dimensionalizing 
            switch (SensorModality)
                case 'RADAR'
                    aux = aux ./ [obj.Re obj.Re/obj.Tc];
                case 'INERTIAL'
                    aux = aux / obj.Re;
            end

            y(1:length(aux),i) = aux;
        else
            y(1,i) = NaN;
        end
    end
end