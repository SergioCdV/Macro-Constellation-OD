
function [Posterior] = CorrectionStep(obj, Measurements, PropPrior, indices)
    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);
    
    % Detection terms 
    L = size(particles,2);
    psi = zeros((length(indices)+1), L);

    % Integration domain 
    nu = obj.dt(1,:);

    for i = 1:length(indices)
        % Extract the observation model and likelihood functions 
        ObservationModel = Measurements{indices(i),5};
        Likelihood = Measurements{indices(i),4};

        for j = 1:L
            % Compute the particle state 
            State = ParticleState(particles(:,j), nu);
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
            psi(i+1,j) = obj.PD * obj.quadrature(l);
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
function [State] = ParticleState(particle, nu)
    % Constants 
    mu = 3.986e14;  % Earth's gravitational parameter
    Re = 6378e3;    % Reference radius of the Earth

    mu = 1; 
    Re = 1; 

    % Delaunay elements of the particle 
    qp = particle(1:4);
    L = particle(5);
    G = particle(6);

    % Solve for the missing Keplerian constants 
    a = Re * L^2;
    e = sqrt(1-(G/L)^2);

    % Preallocation 
    State = zeros(6,length(nu));

    % Compute the radial and velocity distance in the perifocal frame
    r = a*(1-e^2)./(1+e*cos(nu)) .* [cos(nu); sin(nu); zeros(1,length(nu))];
    v = (mu/sqrt(mu * a * (1-e^2))) * [-sin(nu); e + cos(nu); zeros(1,length(nu))];

    Q = QuaternionAlgebra.right_isoclinic( QuaternionAlgebra.quaternion_inverse(qp) );

    for i = 1:length(nu)
        % Compute the mean Brouwer elements from the particle set: Delaunay to mean Brouwer
    
        % Compute the osculating classical orbital elements from the mean Brouwer ones
    
        % Compute the ECI osculating coordinates from the mean Brouwer ones
        aux(:,1) = QuaternionAlgebra.right_isoclinic( [r(:,i); 0] ) * qp;
        aux(:,2) = QuaternionAlgebra.right_isoclinic( [v(:,i); 0] ) * qp;
        aux = Q * aux; 
        State(:,i) = [aux(1:3,1); aux(1:3,2)];
    end

    State = real(State);
end