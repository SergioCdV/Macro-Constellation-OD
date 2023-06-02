

function [planes] = PlanePropagation(obj, planes, step)
    % Constants
    step = step / obj.Tc;
    One = [0;0;0;1];

    % Propagation of the perifocal attitude
    for i = 1:size(planes,2)
        % Action space
        L = planes(5,i);                     % Delaunay action
        G = planes(6,i);                     % Angular momentum
        H = planes(7,i);                     % Polar component of the angular momentum
        eta = G/L;                           % Eccentricity function

        % Compute the angular velocity 
        Omega = zeros(3,1);                  % Angular velocity

        % RAAN motion
        Omega(1) = 3/2*obj.epsilon/(L^7*eta^4) * (H/G);      

        % Exponential mapping
        if (step > obj.Tc / 1e3)
            step = linspace(0,step,1e2);
            dstep = step(2)-step(1);
        else
            dstep = step;
        end

        for j = 1:length(step)
            % Body frame rotation
            Omega = QuaternionAlgebra.right_isoclinic([Omega; 0]) * QuaternionAlgebra.quaternion_inverse(planes(1:4,i)); 
            Omega = QuaternionAlgebra.right_isoclinic(planes(1:4,i)) * Omega; 
            omega = Omega(1:3,1);
    
            % Perigee motion
            omega(3) = omega(3) + 3/4*obj.epsilon/(L^7*eta^4) * (1 - 5*(H/G)^2);            
    
            % Propagate the quaternions only using the Lie-Euler method
            omega = dstep/2 * omega; 
            planes(1:4,i) = QuaternionAlgebra.right_isoclinic(planes(1:4,i)) * QuaternionAlgebra.exp_map([omega; 0], One);
        end
    end
end