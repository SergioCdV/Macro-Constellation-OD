

function [s] = MKV2ECI(mu, x, direction)
    if (direction)
        COE = Astrodynamics.MKV2COE(mu, x, true); 
        s = Astrodynamics.ECI2COE(mu, COE, false).';
    else
        % Compute the angular momentum vector 
        s(1:3,1) = cross(x(1:3,1), x(4:6,1));

        % Compute the eccentricity vector 
        s(4:6,1) = cross(x(4:6), s(1:3,1)) / mu - x(1:3,1) / norm(x(1:3,1));

        % Compute Euler's rotation matrix and the assocaited angles 
        k = s(1:3,1) / norm(s(1:3,1));
        i = s(4:6,1) / norm(s(4:6,1));
        j = cross(k,i); 
        Q = [i j k].';
        Omega = atan2(Q(3,1),-Q(3,2));                      % RAAN
        omega = atan2(Q(1,3),Q(2,3));                       % Argument of perigee

        % Compute the true longitude
        e = norm(s(4:6,1));                                 % Orbital eccentricity 
        r = Q * x(1:3,1);                                   % Position vector in the perifocal frame
        theta = atan2(r(2), r(1));                          % True anomaly
        sin_E = sqrt(1-e^2) * sin(theta)/(1+e*cos(theta));  % Sine of the eccentric anomaly
        cos_E = (cos(theta)+e)/(1+e*cos(theta));            % Cosine of the eccentric anomaly            
        E = atan2(sin_E, cos_E);                            % Eccentric anomaly
        M = E - e * sin(E);                                 % Mean anomaly
        s(7,1) = Omega + omega + M;                         % Longitude
    end
end