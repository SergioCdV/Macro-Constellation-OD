
function [Posterior] = CorrectionStep(obj, Measurements, PropPrior, indices)
    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);
    
    % Corrector step of the KF for the covariance and means

    % Weights update
    L = size(particles,2);
    psi = zeros((length(indices)+1), L);

    for i = 1:length(indices)
        % Extract the observation model and likelihood functions 
        ObservationModel = Measurements{indices(i),5};
        Likelihood = Measurements{indices(i),4};
        SensorModality = Measurements{indices(i),6};

        for j = 1:L
            % Compute the particle state 
            State = ParticleState(obj, particles(:,j));
            l = zeros(1,size(State,2));

            % Evaluate the likelihood
            [~, y] = feval(ObservationModel, State(:,k).');

            if (isempty(y))
                l(k) = 0;
            else
                % Dimensionalizing 
                switch (SensorModality)
                    case 'RADAR'
                        y = y ./ [obj.Re obj.Re/obj.Tc];
                    case 'INERTIAL'
                        y = y / obj.Re;
                end

                l(k) = feval(Likelihood, y.');
            end

            psi(i+1,j) = obj.PD * l;
        end

        psi(i+1,:) = psi(i+1,:) / sum( psi(i+1,:) .* weights(1,:),2 );

        % Check for NaN 
        pos = isnan(psi(i+1,:)); 
        psi(i+1,pos) = zeros(1,sum(pos));
    end

    % Update the weights with the non-detected particles 
    psi(1,:) = (1-obj.PD);
    
    % Assemble the posterior 
    Posterior = [particles; sum(psi,1).*weights];
end

%% Auxiliary functions 
function [State] = ParticleState(obj, particle)
    % Delaunay elements of the particle 
    qp = particle(1:4);
    L = particle(5);
    G = particle(6);

    % Solve for the missing Keplerian constants 
    e = sqrt(1-(G/L)^2);

    % Solve for the true anomaly
    nu = Astrodynamics.KeplerSolver(e,particle(1,1));

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
    State() = [aux(1:3,1); aux(1:3,2)];

    State = real(State);
    State(1:3,:) = obj.Re * State(1:3,:);
    State(4:6,:) = obj.Re/obj.Tc * State(4:6,:);
end