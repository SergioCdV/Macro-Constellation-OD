function [s] = MKV2COE(mu, x, direction)
    if (direction)
        % Compute the Euler angles
        k = x(1:3,1) / norm(x(1:3,1));
        i = x(4:6,1) / norm(x(4:6,1));
        j = cross(k,i); 
        Q = [i j k].';
        Omega = atan2(Q(3,1),-Q(3,2));                      % RAAN
        omega = atan2(Q(1,3),Q(2,3));                       % Argument of perigee
        i = acos(Q(3,3));                                   % Inclination 
        M = x(7,1) - Omega - omega;                         % Mean anomaly
    
        % Compute the geometry of the orbit
        p = dot(x(1:3,1), x(1:3,1), 1) / mu;                % Semilatus rectum
        eta = 1 - dot(x(4:6,1), x(4:6,1));                  % Eccentricity function
        a = p / eta;                                        % Semimajor axis of the orbit 
        e = norm(x(4:6,1));                                 % Orbital eccentricity
    
        % Compute the transformation from COE to ECI
        s = [a e Omega i omega M p];
        s(3:6) = mod(s(3:6),2*pi);
    else
        % Transformation to Cartesian elements 
        s = Astrodynamics.ECI2COE(mu, x, false).';

        % Transformation from ECI to MKV
        s = Astrodynamics.MKV2ECI(mu, s, false);
    end
end