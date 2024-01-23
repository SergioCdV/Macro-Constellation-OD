

function [sc] = ScPropagation(obj, sc, step)
    % Constants
    One = [0;0;0;1];

    if (step < obj.Tc)
        % Propagation of the perifocal attitude
        for i = 1:size(sc,2)
            % Action space
            L = sc(5,i);                     % Delaunay action
            G = sc(6,i);                     % Angular momentum
            H = sc(7,i);                     % Polar component of the angular momentum
            M = sc(8,i);                     % Mean anomaly of the spacecraft
            eta = G/L;                       % Eccentricity function
            cos_i = H/G;                     % Cosine of the inclination
    
            % Compute the angular velocity 
            Omega = zeros(3,1);              % Angular velocity
    
            % RAAN motion
            Omega(1) = 3/2*obj.epsilon/(L^7*eta^4) * (cos_i);     
    
            % Body frame rotation
            Omega = QuaternionAlgebra.RotateVector(sc(1:4,i), Omega); 
            omega = Omega(1:3,1);
    
            % Perigee motion
            omega(3) = omega(3) + 3/4*obj.epsilon/(L^7*eta^4) * (1 - 5 * cos_i^2);            
    
            % Propagate the quaternions only using the Lie-Euler method
            omega = step/2 * omega;
            dq = QuaternionAlgebra.exp_map(omega, One);
            sc(1:4,i) = QuaternionAlgebra.right_isoclinic(dq) * sc(1:4,i);

            % Propagate the anomalies
            omega = 1./L.^3 + (3/4) * obj.epsilon./(L.^4.*eta.^3) .* (1 - 3 * cos_i^2^2);
            sc(8,i) = M + omega * step;
        end
    else
        % Propagation of the perifocal attitude
        for i = 1:size(sc,2)
            % Compute the Delaunay elements
            D = Astrodynamics.Delaunay2MyElements(sc(:,i), false);

            % Action space
            L = D(4,1);                     % Delaunay action
            G = D(5,1);                     % Angular momentum
            H = D(6,1);                     % Polar component of the angular momentum
            eta = G/L;                      % Eccentricity function
            cos_i = H/G;                    % Cosine of the inclination

            % Linear propagation
            D(2,1) = 3/4*obj.epsilon/(L^7*eta^4) * (1 - 5*(cos_i)^2) * step;
            D(3,1) = 3/2*obj.epsilon/(L^7*eta^4) * (cos_i) * step;
            
            % Re-assembling
            aux = Astrodynamics.Delaunay2MyElements(D, true);
            sc(1:4,i) = aux(1:4);
        end
    end
end