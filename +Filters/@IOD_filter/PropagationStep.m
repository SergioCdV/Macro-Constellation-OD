
function [PropPrior] = PropagationStep(obj, last_epoch, new_epoch, Prior)
    % Constants 
    epsilon = -1.08263e-3;              % J2
    One = [0;0;0;1];                    % Identity quaternion
    mu = 3.986e14;                      % Earth's gravitational parameter
    Re = 6378e3;                        % Reference radius of the Earth
    Tc = sqrt(Re^3/mu);                 % Characteristic time

    step = new_epoch - last_epoch;      % Time step
    step = step / Tc;                   % Dimensionalizing

    % Preallocation 
    particles = Prior(1:end-1,:);       % Prior states
    weights = Prior(end,:);             % Prior weights
    
    % Perform the propagation 
    if (true && step > 0)
        for i = 1:size(particles,2)
            % Extract the actions 
            L = particles(5,i);                     % Delaunay action
            G = particles(6,i);                     % Angular momentum
            H = particles(7,i);                     % Polar component of the angular momentum
            eta = G/L;                              % Eccentricity function

            % Compute the angular velocity 
            Omega = zeros(3,1);                 % Angular velocity

            % RAAN motion
            Omega(1) = 3/2*epsilon/(L^7*eta^2) * (H/G);                        

            % Body frame rotation
            Omega = QuaternionAlgebra.right_isoclinic([Omega; 0]) * QuaternionAlgebra.quaternion_inverse(particles(1:4,i)); 
            Omega = QuaternionAlgebra.right_isoclinic(particles(1:4,i)) * Omega; 
            omega = Omega(1:3,1);

            % Perigee motion
            omega(3) = omega(3) +  3/4*epsilon/(L^7*eta^3) * (1 - 5*(H/G)^2);            

            % Propgate the quaternions only using the Lie-Euler method
            omega = step/2 * omega; 
            particles(1:4,i) = QuaternionAlgebra.right_isoclinic(particles(1:4,i)) * QuaternionAlgebra.exp_map([omega; 0], One);
        end
    end

    % Assemble the final propagated prior 
    PropPrior = [particles; weights];
end