
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
    psi = zeros(1,(length(indices)+1) * L);

    particles = repmat(particles, 1, length(indices)+1); 
    weights = repmat(weights, 1, length(indices)+1); 

    for i = 1:length(indices)
        % Extract the observation model and likelihood functions 
        Z = Measurements{indices(i),2}(2:end).';
        Likelihood = Measurements{indices(i),4};
        ObservationModel = Measurements{indices(i),5};
        R = reshape(Measurements{indices(i),6}, [size(Z,1) size(Z,1)]);
        SensorModality = Measurements{indices(i),7};
        REstimator.R = R;
        REstimator = REstimator.AssignObservationProcess(size(Z,1), @(state)FullProcess(obj, ObservationModel, state));

        for j = 1:L
            index = 1 + L * (i-1) + j;
            sigma_points = reshape( particles(pos^2+pos+1:end,j), pos, []);
            Sigma = reshape( particles(pos+1:pos^2+pos,j), [pos pos]);
            [mu, S, ~, y] = REstimator.CorrectionStep(sigma_points, particles(1:pos,j), Sigma, Z);
            particles(1:pos,index) = mu;
            particles(pos+1:pos+pos^2,index) = reshape(S, [], 1);
            mu

            if (isempty(y))
                l = 0;
            else
                % Dimensionalizing 
                switch (SensorModality)
                    case 'RADAR'
                        y = y ./ [obj.Re; obj.Re/obj.Tc];
                    case 'INERTIAL'
                        y = y / obj.Re;
                end

                l = feval(Likelihood, y);
            end

            psi(1, 1 + L*(i-1)+ j) = obj.PD * l;
        end

        psi(1,1+L*(i-1):L*i) = psi(1,1+L*(i-1):L*i) / sum( psi(1,1+L*(i-1):L*i) .* weights(1,1+L*(i-1):L*i),2 );

        % Check for NaN 
        pos_nan = isnan(psi(1,1+L*(i-1):L*i)); 
        psi(1,pos_nan) = zeros(1,sum(pos_nan));
    end

    % Update the weights with the non-detected particles 
    psi(1,1:L) = (1-obj.PD);
    
    % Assemble the posterior 
    Posterior = [particles([1:pos+pos^2],:); dot(psi,weights,1)];
end

%% Auxiliary functions 
% Full measurement process from the Delaunay set 
function [y] = FullProcess(obj, ObservationModel, state)
    % Preallocation 
    for i = 1:size(state,2)
        State = ParticleState(obj, state(:,i)).';
        [~, y(:,i)] = feval(ObservationModel, State);
    end
end

% Transformation from Delaunay to Cartesian elements
function [State] = ParticleState(obj, particle)
    % Delaunay elements of the particle 
    qp = particle(1:4,1);
    L = particle(5,1);
    G = particle(6,1);

    % Solve for the missing Keplerian constants 
    e = sqrt(1-(G/L)^2);

    % Solve for the true anomaly
    nu = Astrodynamics.KeplerSolver(e,particle(7,1));

    % Compute the radial and velocity distance in the perifocal frame
    r = G./(1+e*cos(nu)) .* [cos(nu); sin(nu); zeros(1,length(nu))];
    v = [-sin(nu); e + cos(nu); zeros(1,length(nu))] / G;

    Q = QuaternionAlgebra.right_isoclinic( QuaternionAlgebra.quaternion_inverse(qp) );

    % Compute the mean Brouwer elements from the particle set: Delaunay to mean Brouwer

    % Compute the osculating classical orbital elements from the mean Brouwer ones

    % Compute the ECI osculating coordinates from the mean Brouwer ones
    aux(:,1) = QuaternionAlgebra.right_isoclinic( [r; 0] ) * qp;
    aux(:,2) = QuaternionAlgebra.right_isoclinic( [v; 0] ) * qp;
    aux = Q * aux; 
    State = [aux(1:3,1); aux(1:3,2)];

    State = real(State);
    State(1:3,:) = obj.Re * State(1:3,:);
    State(4:6,:) = obj.Re/obj.Tc * State(4:6,:);
end