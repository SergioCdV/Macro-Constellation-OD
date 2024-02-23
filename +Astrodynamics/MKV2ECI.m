function [s] = MKV2ECI(mu, x, direction)
    if (direction)
        COE = Astrodynamics.MKV2COE(mu, x, true); 
        s = Astrodynamics.ECI2COE(mu, COE, false);
    else
        % Preallocation 
        s = zeros(6,size(x,2));     % Transformed state

        % Compute the angular momentum vector 
        s(1:3,:) = cross(x(1:3,:), x(4:6,:));

        % Compute the eccentricity vector 
        s(4:6,:) = cross(x(4:6,:), s(1:3,:)) / mu - x(1:3,:) ./ sqrt( dot(x(1:3,:), x(1:3,:), 1) );

        % Compute Euler's rotation matrix and the assocaited angles 
        k = s(1:3,:) ./ sqrt( dot(s(1:3,:), s(1:3,:), 1) );
        I = s(4:6,:) ./ sqrt( dot(s(4:6,:), s(4:6,:), 1) );
        j = cross(k,I); 

        % Transformation
        Omega = zeros(1,size(x,2)); 
        omega = Omega;
        theta = Omega;

        for i = 1:size(x,2)
            Q = [I(:,i) j(:,i) k(:,i)].';
            Omega(i) = atan2(Q(3,1),-Q(3,2));                         % RAAN
            omega(i) = atan2(Q(1,3),Q(2,3));                          % Argument of perigee
    
            % Compute the true longitude
            r = Q * x(1:3,i);                                         % Position vector in the perifocal frame
            theta(i) = atan2(r(2), r(1));                             % True anomaly
        end

        e = sqrt( dot(s(4:6,:), s(4:6,:), 1) );                       % Orbital eccentricity 
        sin_E = sqrt(1-e.^2) .* sin(theta) ./ (1 + e .* cos(theta));  % Sine of the eccentric anomaly
        cos_E = (cos(theta) + e) ./ (1 + e .* cos(theta));            % Cosine of the eccentric anomaly            
        E = atan2(sin_E, cos_E);                                      % Eccentric anomaly
        M = E - e .* sin(E);                                          % Mean anomaly
        s(7,:) = Omega + omega + M;                                   % Longitude
    end
end