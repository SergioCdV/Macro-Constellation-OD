
function [Posterior] = CorrectionStep(obj, Likelihood, ObservationModel, PropPrior)
    % Extract elements 
    particles = PropPrior(1:end-1,:);
    weights = PropPrior(end,:);
    
    % Detection terms 
    psi = zeros(size(Likelihood,1),size(weights,2));
    G = psi;
    for i = 1:size(Likelihood,1)
        for j = 1:size(weights,2)
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
            psi(i,j) = obj.PD * max(l);
        end
        G(i,:) = psi(i,:) / sum( psi(i,:) .* weights );
    end

    % Update the weights with the likelihood function 
    for i = 1:size(weights,2)
        % Misdetection term 
        Kmis = 1-obj.PD;   

        % Detection term
        Kdet = sum(G,1);
        index = isnan(Kdet);
        if (any(index))
            Kdet(index) = zeros(1,length(index));
        end

        % Update
        K = Kmis + Kdet; 
        weights(i) = K(i) * weights(i);
    end

    % Assemble the posterior 
    Posterior = [particles; weights];
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