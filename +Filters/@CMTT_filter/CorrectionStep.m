
function [Posterior] = CorrectionStep(obj, indices, Measurements, Estimator, PropPrior)
    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);

    pos = Estimator.StateDim+1;
    Estimator.StateDim = 4; 
    REstimator = Estimator.Init();
    
    % Corrector step of the KF for the covariance and means
    Sigma = reshape(particles(pos+1:pos+(pos-1)^2,1), [pos-1 pos-1]); 
    F = Sigma(1:3,1:3)^(-1);

    % Weights update
    L = size(particles,2);
    psi = zeros(1, L * (length(indices) + 1));
    particles = repmat(particles, 1, length(indices)+1); 

    for i = 1:length(indices)
        % Extract the observation model and likelihood functions 
        Z = Measurements{indices(i),2}(2:end).';
        Likelihood = Measurements{indices(i),4};
        ObservationModel = Measurements{indices(i),5};
        R = reshape(Measurements{indices(i),6}, [size(Z,1) size(Z,1)]);
        SensorModality = Measurements{indices(i),7};
        REstimator.R = R;
        REstimator = REstimator.AssignObservationProcess(size(Z,1), @(state)FullProcess(obj, SensorModality, ObservationModel, state));

        for j = 1:L
            index = L * i + j;
            Sigma = blkdiag(1, reshape( particles(pos+1:pos+(pos-1)^2,j), [pos-1 pos-1]));
            sigma_points = reshape( particles((pos-1)^2+pos+1:end,j), pos, []);

            % UKF step
            [mu, S, ~, y] = REstimator.CorrectionStep(sigma_points, particles(1:pos,j), Sigma, Z);
            Sigma = Sigma(2:end,2:end);
            S = S(2:end,2:end);

            % Deconditioning
            mu(5:end) = mu(5:end) - Sigma(4:end,1:3) * F * obj.Gibbs_vector(:,i);
            S(4:end,4:end) = S(4:end,4:end) + S(4:end,1:3) * F * S(4:end,1:3).';
            S = 0.5 * (S + S.') + obj.PD_tol * eye(size(S,1));

            % Particles update
            particles(1:pos,index) = mu;
            particles(pos+1:pos+(pos-1)^2,index) = reshape(S, [], 1);

            if (isempty(y) || any(isnan(y)))
                l = 0;
            else
                l = feval(Likelihood, y);
            end

            psi(1, L*i + j) = obj.PD * l;
        end

        psi(1,1+L*i:L*(i+1)) = psi(1,1+L*i:L*(i+1)) .* weights(1,:) / sum( psi(1,1+L*i:L*(i+1)) .* weights(1,:),2 );
    end

    % Update the weights with the non-detected particles 
    psi(1,1:L) = (1-obj.PD) .* weights(1,:);

    % Check for NaN 
    pos_nan = isnan(psi(1,:)); 
    psi(1,pos_nan) = zeros(1,sum(pos_nan));
    
    % Assemble the posterior 
    Posterior = [particles(1:pos+(pos-1)^2,:); psi];
end

%% Auxiliary functions 
% Full measurement process from the Delaunay set 
function [y] = FullProcess(obj, SensorModality, ObservationModel, state)
    % Preallocation 
    for i = 1:size(state,2)
        State = obj.ParticleState(SensorModality, state(:,i)).';
        [~, aux] = feval(ObservationModel, State);
        if (~isempty(aux))
            % Dimensionalizing 
            switch (SensorModality)
                case 'RADAR'
                    aux = aux ./ [obj.Re obj.Re/obj.Tc];
                case 'INERTIAL'
                    aux = aux / obj.Re;
            end

            y(1:length(aux),i) = aux.';
        else
            y(1,i) = NaN;
        end
    end
end