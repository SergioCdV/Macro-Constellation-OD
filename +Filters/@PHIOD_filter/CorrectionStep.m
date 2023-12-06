
function [Posterior, mean_PD] = CorrectionStep(obj, indices, Measurements, Estimator, PropPrior)
    % Constants 
    pos = 8; 

    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);
    
    % Weights update
    L = size(particles,2);
    psi = zeros(1, L * (length(indices) + 1));
    particles = repmat(particles, 1, length(indices)+1); 
    PD_t = zeros(1,L);

    % Box constraints 
    box(1,:) = [obj.Lmin obj.Lmax];
    box(2,:) = [sqrt(1-obj.emax^2) 1];
    box(3,:) = [-1 1];
       
    for i = 1:length(indices)
        % Extract the observation model and likelihood functions 
        Z = Measurements{indices(i),2}(2:end).';
        Likelihood = Measurements{indices(i),4};
        ObservationModel = Measurements{indices(i),5};
        R = reshape(Measurements{indices(i),6}, [size(Z,1) size(Z,1)]);
        SensorModality = Measurements{indices(i),7};
        Estimator.R = R;

        Estimator = Estimator.AssignObservationProcess(size(Z,1), @(state)FullProcess(obj, SensorModality, ObservationModel, box, state));

        for j = 1:L
            % Updated particle index
            index = L * i + j;

            % UKF step
            mu = particles(1:pos+3,j);                                % Predicted mean
            Sigma = particles(pos+3+1:pos+3+(pos-1)^2,j);             % Predicted covariance
            Sigma = reshape(Sigma, [pos-1 pos-1]);
            sigma_points = particles(pos+3+(pos-1)^2+1:end,j);        % Predicted sigma points
            sigma_points = reshape(sigma_points, pos+3, []);
        
            [mu, S, ~, y] = Estimator.CorrectionStep(sigma_points, mu, Sigma, Z);

            if (isempty(y) || any(isnan(y)))
                l = 0;
                PD_t(j) = 0;
            else
                % ADMM step 
                X = projection_constraints(mu(4:6,1), box);
                R = QuaternionAlgebra.Quat2Matrix(mu(pos:pos+3,:));
                X(end) = X(end-1) * R(3,3);

                % Final assembling
                particles(1:pos+3,index) = [mu(1:3,1); X; mu(pos-1:pos+3,:)];
                particles(pos+3+1:pos+3+(pos-1)^2,index) = reshape(S, [], 1);

                % Likelihood
                l = feval(Likelihood, y);
                PD_t(j) = obj.PD;
            end
            
            psi(1, L*i + j) = PD_t(j) * l;
        end

        psi(1,1+L*i:L*(i+1)) = psi(1,1+L*i:L*(i+1)) .* weights(1,:) / sum( psi(1,1+L*i:L*(i+1)) .* weights(1,:),2 );
    end

    % Update the weights with the non-detected particles 
    psi(1,1:L) = (1-PD_t) .* weights(1,:);

    % Check for NaN 
    pos_nan = isnan(psi(1,:)); 
    psi(1,pos_nan) = zeros(1,sum(pos_nan));
    
    % Assemble the posterior 
    particles = [particles(pos:pos+3,:); particles(4:pos-1,:); particles(pos+3+1:pos+3+(pos-1)^2,:); ];
    Posterior = [particles; psi];

    mean_PD = mean(PD_t);
end

%% Auxiliary functions 
% Full measurement process from the Delaunay set 
function [y] = FullProcess(obj, SensorModality, ObservationModel, box, particles)
    % Compute the spacecraft state
    particles(5:7,:) = projection_constraints(particles(5:7,:), box);
    State = real(obj.ParticleState(SensorModality, particles));

    for i = 1:size(State,2)
        [~, aux] = feval(ObservationModel, State(:,i).');
        if (~isempty(aux) && all(~isnan(aux)))
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

% Projection function 
function [X] = projection_constraints(X, box)
    for k = 1:size(X,2)
        % Projection of the Delaunay action 
        if (X(1,k) < box(1,1))
            X(1,k) = box(1,1);
        elseif (X(1,k) > box(1,2))
            X(1,k) = box(1,2);
        end

        if (X(1,k) == 0)
            a = 1;
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