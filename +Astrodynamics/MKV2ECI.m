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
        % Transformation to COE
        COE = Astrodynamics.MKV2COE(mu, x, true); 

        % Transformation to ECI
        s = Astrodynamics.ECI2COE(mu, COE, false);
    else
        % Constants 
        r = x(1:3,:);               % Inertial position vector 
        v = x(4:6,:);               % Inertial velocity vector 

        % Preallocation 
        s = zeros(6,size(x,2));     % Transformed state

        % Compute the angular momentum vector 
        s(1:3,:) = cross(r, v);

        % Compute the eccentricity vector 
        s(4:6,:) = cross(v, s(1:3,:)) / mu - r ./ sqrt( dot(r, r, 1) );     % Eccentricity vector
        e = sqrt( dot(s(4:6,:), s(4:6,:), 1) );                             % Orbital eccentricity 

        % Compute Euler's rotation matrix and the associated angles 
        k = s(1:3,:) ./ sqrt( dot(s(1:3,:), s(1:3,:), 1) );
        
        I = s(4:6,:) ./ e;
        K = repmat([0;0;1], 1, size(x,2));                     % Inertial Z axis
        n = cross(K, s(1:3,:));                                % Node vector
        I(:, e == 0) = n(:, e == 0);                           % Circular orbit singularity

        j = cross(k,I); 

        % Position in the perifocal frame 
        Q = reshape([I; j; k], 3, []).';
        Omega = atan2(Q(3:3:end,1),-Q(3:3:end,2));             % RAAN
        omega = atan2(Q(1:3:end,3),Q(2:3:end,3));              % Argument of perigee
        i = acos(Q(3:3:end,3));                                % Inclination

        r0 = r;
        for j = 1:size(r,2)
            r0(:,j) = Q(1+3*(j-1):3*j,:) * r(:,j);
        end

        theta = atan2(r0(2,:), r0(1,:));                              % True anomaly

        % Longitude
        den = 1 + e .* cos(theta);
        sin_E = sqrt(1-e.^2) .* sin(theta) ./ den;                    % Sine of the eccentric anomaly
        cos_E = (cos(theta) + e) ./ den;                              % Cosine of the eccentric anomaly            
        E = atan2(sin_E, cos_E);                                      % Eccentric anomaly
        M = E - e .* sin(E);                                          % Mean anomaly

        % Non-singular COE
        elements = [zeros(1,size(e,2)); e; Omega.'; i.'; omega.'; M];
        elements = Astrodynamics.rv_singularity(s(4:6,:), n, r, Q, elements);  
        
        s(7,:) = Omega.' + omega.' + elements(6,:);                   % Mean longitude
    end
end