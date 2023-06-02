
function [Posterior] = CorrectionStep(obj, indices, Measurements, Estimator, PropPrior)
    % Constants 
    pos = 7; 

    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);
    
    % Weights update
    L = size(particles,2);
    psi = zeros(1, L * (length(indices) + 1));
    particles = repmat(particles, 1, length(indices)+1); 

    % Optimization options 
    options = optimoptions('fmincon', 'Display', 'none');

    for i = 1:length(indices)
        % Extract the observation model and likelihood functions 
        Z = Measurements{indices(i),2}(2:end).';
        Likelihood = Measurements{indices(i),4};
        ObservationModel = Measurements{indices(i),5};
        R = reshape(Measurements{indices(i),6}, [size(Z,1) size(Z,1)]);
        SensorModality = Measurements{indices(i),7};
        Estimator.R = R;

        for j = 1:L
            index = L * i + j;
            Sigma = reshape( particles(pos+3+1:pos+3+(pos-1)^2,j), [pos-1 pos-1]);
            sigma_points = reshape( particles(pos+3+(pos-1)^2+1:end,j), pos+3, []);

            % Solve for the ML anomaly
            mu = [particles(pos:pos+3,j); particles(4:6,j)];
            theta = fmincon(@(theta)LikeProcess(obj, SensorModality, ObservationModel, Likelihood, mu, theta), 0, [], [], [], [], 0, 2*pi, [], options); 

            % UKF step
            Estimator = Estimator.AssignObservationProcess(size(Z,1), @(state)FullProcess(obj, SensorModality, ObservationModel, state, theta));
            [mu, S, ~, ~] = Estimator.CorrectionStep(sigma_points, particles(1:pos+3,j), Sigma, Z);

            % Particles update
            particles(1:pos+3,index) = mu;
            particles(pos+3+1:pos+3+(pos-1)^2,index) = reshape(S, [], 1);

            % Weights update
            fun = @(theta)LikeProcess(obj, SensorModality, ObservationModel, Likelihood, [mu(7:end,1); mu(4:6,1)], theta);
            l = integral(fun, 0, 2*pi);
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
    particles = [particles(pos:pos+3,:); particles(4:pos-1,:); particles(pos+3+1:pos+3+(pos-1)^2,:); ];
    Posterior = [particles; psi];
end

%% Auxiliary functions 
% Full measurement process from the Delaunay set 
function [y] = FullProcess(obj, SensorModality, ObservationModel, particles, M)
    % Compute the spacecraft state
    State = obj.ParticleState(SensorModality, particles, M);

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

% Full measurement process from the Delaunay set 
function [l] = LikeProcess(obj, SensorModality, ObservationModel, Likelihood, particles, M)
    % Preallocation 
    l = zeros(1,length(M)); 

    % Valuation of the likelihood
    for i = 1:length(M)
        State = obj.ParticleState(SensorModality, particles, M(i));
        [~, aux] = feval(ObservationModel, State.');
        if (~isempty(aux) && all(~isnan(aux)))
            % Dimensionalizing 
            switch (SensorModality)
                case 'RADAR'
                    aux = aux ./ [obj.Re obj.Re/obj.Tc];
                case 'INERTIAL'
                    aux = aux / obj.Re;
            end

            l(i) = feval(Likelihood, aux.');
        else
            l(i) = 0;
        end
    end
end