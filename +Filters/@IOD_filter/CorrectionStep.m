
function [Posterior] = CorrectionStep(obj, Measurements, PropPrior, indices)
    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);
    
    % Detection terms 
    L = size(particles,2);
    psi = zeros(1, (length(indices)+1) * L);

    for i = 1:length(indices)
        % Extract the observation model and likelihood functions 
        ObservationModel = Measurements{indices(i),5};
        Likelihood = Measurements{indices(i),4};

        for j = 1:L
            % Compute the particle state 
            State = ParticleState(particles(:,j));
            l = zeros(1,size(State,2));

            for k = 1:size(State,2)
                % Evaluate the likelihood
                [~, y] = feval(ObservationModel, State(:,k).');
                if (isempty(y))
                    l(k) = 0;
                else
                    l(k) = feval(Likelihood, y.');
                end
            end
            psi(1, i * L + j) = obj.PD * max(l) * weights(1,j);
        end
        psi(1, i*L:(i+1)*L) = psi(1,i*L:(i+1)*L) / sum( psi(1,i*L:(i+1)*L) );
    end

    % Check for NaN 
    pos = isnan(psi); 
    psi(1,pos) = zeros(1,sum(pos));

    % Update the weights with the non-detected particles 
    psi(1,1:L) = (1-obj.PD) * weights;
    
    % Assemble the posterior 
    Posterior = [repmat(particles, 1, length(indices)+1); psi];
end

%% Auxiliary functions 
function [State] = ParticleState(particle)
    % Constants 
    mu = 3.986e14;  % Earth's gravitational parameter
    Re = 6378e3;    % Reference radius of the Earth

    % Delaunay elements of the particle 
    qp = particle(1:4);
    L = particle(5);
    G = particle(6);

    % True anomaly space
    nu = linspace(0,2*pi,1e2);

    % Solve for the missing Keplerian constants 
    a = Re * L^2;
    e = sqrt(1-(G/L)^2);

    % Preallocation 
    State = zeros(6,length(nu));

    % Compute the radial and velocity distance in the perifocal frame
    r = a*(1-e^2)./(1+e*cos(nu)) .* [cos(nu); sin(nu); zeros(1,length(nu))];
    v = (mu/sqrt(mu * a * (1-e^2))) * [-sin(nu); e + cos(nu); zeros(1,length(nu))];

    for i = 1:length(nu)
        % Compute the mean Brouwer elements from the particle set: Delaunay to mean Brouwer
    
        % Compute the osculating classical orbital elements from the mean Brouwer ones
    
        % Compute the ECI osculating coordinates from the mean Brouwer ones
        aux = QuaternionAlgebra.right_isoclinic( QuaternionAlgebra.quaternion_inverse(qp) ) * (QuaternionAlgebra.right_isoclinic( [r(:,i); 0] ) * qp);
        State(1:3,i) = aux(1:3); 
        aux = QuaternionAlgebra.right_isoclinic( QuaternionAlgebra.quaternion_inverse(qp) ) * (QuaternionAlgebra.right_isoclinic( [v(:,i); 0] ) * qp);
        State(4:6,i) = aux(1:3); 
    end

    State = real(State);
end