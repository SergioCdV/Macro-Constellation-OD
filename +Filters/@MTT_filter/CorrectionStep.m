
function [Posterior] = CorrectionStep(obj, indices, Measurements, Estimator, PropPrior)
    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);

    pos = Estimator.StateDim;
    Estimator.StateDim = 4; 
    REstimator = Estimator.Init();
    
    % Corrector step of the KF for the covariance and means

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
            sigma_points = reshape( particles(pos^2+pos+1:end,j), pos, []);
            Sigma = reshape( particles(pos+1:pos^2+pos,j), [pos pos]);
            [mu, S, ~, y] = REstimator.CorrectionStep(sigma_points, particles(1:pos,j), Sigma, Z);

            if (any(isnan(S)))
                a = 1;
            end
            particles(1:pos,index) = mu;
            particles(pos+1:pos+pos^2,index) = reshape(S, [], 1);

            if (isempty(y) || any(isnan(y)))
                l = 0;
            else
                l = feval(Likelihood, y);
            end

            psi(1, L*i + j) = obj.PD * l;
        end

        psi(1,1+L*i:L*(i+1)) = psi(1,1+L*i:L*(i+1)) .* weights(1,:) / sum( psi(1,1+L*i:L*(i+1)) .* weights(1,:),2 );
    end

    % Check for NaN 
    pos_nan = isnan(psi(1,:)); 
    psi(1,pos_nan) = zeros(1,sum(pos_nan));

    % Update the weights with the non-detected particles 
    psi(1,1:L) = (1-obj.PD) .* weights(1,:);
    
    % Assemble the posterior 
    Posterior = [particles(1:pos+pos^2,:); psi];
end

%% Auxiliary functions 
% Full measurement process from the Delaunay set 
function [y] = FullProcess(obj, SensorModality, ObservationModel, state)
    % Preallocation 
    for i = 1:size(state,2)
        State = ParticleState(obj, state(:,i)).';
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

% Transformation from Delaunay to Cartesian elements
function [State] = ParticleState(obj, particle)
    % Delaunay elements of the particle 
    qp = particle(1:4,1);
    L = particle(5,1);
    G = particle(6,1);
    H = particle(7,1);
    M = particle(8,1);
 
    % Compute the RAAN and AoP from qp 
    diff = atan2(qp(2,1), qp(1,1));
    plus = atan2(qp(3,1), qp(4,1));
    Omega = 2 * (plus+diff);
    omega = 2 * plus - Omega;

    % Assemble the set
    D = [M Omega omega L G H].';     

    % Compute the ECI osculating coordinates from the long-period ones
    % Do = Astrodynamics.Brouwer_solution(obj.epsilon, D);
%     Do = D;
% 
%     % Preallocation 
%     State = zeros(6,size(D,2)); 
%     for i = 1:size(D,2)
%         State(:,i) = Astrodynamics.Delaunay2ECI(Do);
%         State(1:3,i) = obj.Re * State(1:3,i);
%         State(4:6,i) = obj.Re/obj.Tc * State(4:6,i);
%     end

    State = zeros(6,size(D,2)); 
    for i = 1:size(D,2)
        State(:,i) = Astrodynamics.Lara_solution(obj.epsilon, D);
        State(1:3,i) = obj.Re * State(1:3,i);
        State(4:6,i) = obj.Re/obj.Tc * State(4:6,i);
    end
end