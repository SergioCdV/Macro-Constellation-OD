%% Constellation tracking
% Author: Sergio Cuevas
% Date: 13/03/2024

%% Milankovitch's to ECI transformation
% The following function provides the transformation between Milankovitch's
% element set (angular momentum, LRL vector and true longitude) and the
% Cartesian position and velocity state

% Inputs: -mu, scalar, the gravitational constant of the central body
%         -x, an array of 6 x n, containing the set of states to be transformed along the column space
%         -direction, a boolean to indicate the direction of the transformation (1 for MKV 2 ECI, else for ECI 2 MKV)

% Ouput: - s, an array with the same dimensions of x, containing the transformed variables

function [s] = MKV2ECI(mu, x, direction)
    if (direction)
        COE = utils.MKV2COE(mu, x, true); 
        s = utils.ECI2COE(mu, COE, false);
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