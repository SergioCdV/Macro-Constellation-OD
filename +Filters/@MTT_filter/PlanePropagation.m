

function [planes] = PlanePropagation(obj, last_epoch, prop_epoch, planes)
    % Constants
    step = prop_epoch - last_epoch;         % Time step 
    step = 86400 / obj.Tc * step;
    One = [0;0;0;1];

    if (step < obj.Tc)
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
    
            % Body frame rotation
            Omega = QuaternionAlgebra.right_isoclinic([Omega; 0]) * QuaternionAlgebra.quaternion_inverse(planes(1:4,i)); 
            Omega = QuaternionAlgebra.right_isoclinic(planes(1:4,i)) * Omega; 
            omega = Omega(1:3,1);
    
            % Perigee motion
            omega(3) = omega(3) + 3/4*obj.epsilon/(L^7*eta^4) * (1 - 5*(H/G)^2);            
    
            % Propagate the quaternions only using the Lie-Euler method
            omega = step/2 * omega; 
            dq = QuaternionAlgebra.exp_map([omega; 0], One);
            planes(1:4,i) = QuaternionAlgebra.right_isoclinic(dq) * planes(1:4,i);
        end
    else
        % Propagation of the perifocal attitude
        for i = 1:size(planes,2)
            % Action space
            L = planes(5,i);                     % Delaunay action
            G = planes(6,i);                     % Angular momentum
            H = planes(7,i);                     % Polar component of the angular momentum
            eta = G/L;                           % Eccentricity function

            % Compute the Delaunay elements
            D = Astrodynamics.Delaunay2MyElements(planes(:,i), false);

            % Linear propagation
            D(2,1) = 3/4*obj.epsilon/(L^7*eta^4) * (1 - 5*(H/G)^2) * step;
            D(3,1) = 3/2*obj.epsilon/(L^7*eta^4) * (H/G) * step;
            
            % Re-assembling
            aux = Astrodynamics.Delaunay2MyElements(D, true);
            planes(1:4,i) = aux(1:4);
        end
    end
end