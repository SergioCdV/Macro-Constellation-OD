
function [PropPrior] = PropagationStep(obj, last_epoch, new_epoch, Prior)
    % Constants 
    One = [0;0;0;1];                    % Identity quaternion
    step = new_epoch - last_epoch;      % Time step
    step = 86400 / obj.Tc * step;

    if (step > obj.Tc / 1e3)
        step = linspace(0,step,1e2);
        dstep = step(2)-step(1);
    else
        dstep = step;
    end

    % Preallocation 
    particles = Prior(1:end-1,:);       % Prior states
    weights = Prior(end,:);             % Prior weights
    
    % Perform the propagation 
    for i = 1:size(particles,2)
        % Action space
        L = particles(5,i);                  % Delaunay action
        G = particles(6,i);                  % Angular momentum
        H = particles(7,i);                  % Polar component of the angular momentum
        eta = G/L;                           % Eccentricity function

        % Compute the angular velocity 
        Omega = zeros(3,1);                  % Angular velocity

        % RAAN motion
        Omega(1) = 3/2*obj.epsilon/(L^7*eta^4) * (H/G);      

        % Exponential mapping
        for j = 1:length(step)
            % Body frame rotation
            Omega = QuaternionAlgebra.right_isoclinic([Omega; 0]) * QuaternionAlgebra.quaternion_inverse(particles(1:4,i)); 
            Omega = QuaternionAlgebra.right_isoclinic(particles(1:4,i)) * Omega; 
            omega = Omega(1:3,1);
    
            % Perigee motion
            omega(3) = omega(3) +  3/4*obj.epsilon/(L^7*eta^4) * (1 - 5*(H/G)^2);            
    
            % Propgate the quaternions only using the Lie-Euler method
            omega = dstep/2 * omega; 
            dq = QuaternionAlgebra.exp_map([omega; 0], One);
            planes(1:4,i) = QuaternionAlgebra.right_isoclinic(dq) * planes(1:4,i);
        end
    end

    % Assemble the final propagated prior 
    PropPrior = [particles; weights];
end