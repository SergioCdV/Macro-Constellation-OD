function [s] = MKV2COE(mu, x, direction)
    
    if (direction)
        k = x(1:3,:) ./ sqrt( dot(x(1:3,:), x(1:3,:), 1) );
        I = x(4:6,:) ./ sqrt( dot(x(4:6,:), x(4:6,:), 1) );
        j = cross(k,I); 

        for i = 1:size(x,2)
            % Compute the Euler angles
            Q = [I(:,i) j(:,i) k(:,i)].';
            Omega = atan2(Q(3,1),-Q(3,2));                  % RAAN
            omega = atan2(Q(1,3),Q(2,3));                   % Argument of perigee
            I = acos(Q(3,3));                               % Inclination 
        end

        M = x(7,:) - Omega - omega;                         % Mean anomaly
    
        % Compute the geometry of the orbit
        p = dot(x(1:3,:), x(1:3,:), 1) / mu;                % Semilatus rectum
        eta = 1 - dot(x(4:6,:), x(4:6,:), 1);               % Eccentricity function
        a = p ./ eta;                                       % Semimajor axis of the orbit 
        e = sqrt( dot(x(4:6,:), x(4:6,:), 1) );             % Orbital eccentricity
    
        % Compute the transformation from COE to ECI
        s = [a; e; Omega; I; omega; M; p];
        s(3:6,:) = mod(s(3:6,:), 2 * pi);
    else
        % Transformation to Cartesian elements 
        s = Astrodynamics.ECI2COE(mu, x, false);

        % Transformation from ECI to MKV
        s = Astrodynamics.MKV2ECI(mu, s, false);
    end
end