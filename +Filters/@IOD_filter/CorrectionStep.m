
function [Posterior] = CorrectionStep(obj, indices, Measurements, PropPrior)
    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);
    
    % Detection terms 
    L = size(particles,2);
    psi = zeros((length(indices)+1), L);
    l = zeros(1,size(obj.nu,2));

    for i = 1:length(indices)
        % Extract the observation model and likelihood functions 
        ObservationModel = Measurements{indices(i),5};
        Likelihood = Measurements{indices(i),4};
        SensorModality = Measurements{indices(i),7};

        for j = 1:L
            % Compute the particle state 
            State = obj.ParticleState(SensorModality, particles(:,j), obj.nu);
            [y] = FullProcess(obj, SensorModality, ObservationModel, State);

            % Evaluate the likelihood
            for k = 1:size(obj.nu, 2)
                if (isempty(y(:,k)) || any(isnan(y(:,k))))
                    l(k) = 0;
                else
                    l(k) = feval(Likelihood, y(:,k));
                end
            end

            psi(i+1,j) = obj.PD * obj.quadrature(l);
        end

        psi(i+1,:) = psi(i+1,:) / sum( psi(i+1,:) .* weights(1,:), 2);

        % Check for NaN 
        pos = isnan(psi(i+1,:)); 
        psi(i+1,pos) = zeros(1,sum(pos));
    end

    % Update the weights with the non-detected particles 
    psi(1,:) = (1-obj.PD);

    % Check for NaN 
    pos = isnan(psi(1,:)); 
    psi(1,pos) = zeros(1,sum(pos));
    
    % Assemble the posterior 
    Posterior = [particles; sum(psi,1).*weights];
end

%% Auxiliary functions 
% Full measurement process from the Delaunay set 
function [y] = FullProcess(obj, SensorModality, ObservationModel, State)
    % Preallocation 
    for i = 1:size(State,2)
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