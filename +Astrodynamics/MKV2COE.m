%% Constellation tracking
% Author: Sergio Cuevas
% Date: 13/03/2024

%% Milankovitch's to COE transformation
% The following function provides the transformation between Milankovitch's
% element set (angular momentum, LRL vector and true longitude) and the
% classical orbital elements

% Inputs: -mu, scalar, the gravitational constant of the central body
%         -x, an array of 6 x n, containing the set of states to be transformed along the column space
%         -direction, a boolean to indicate the direction of the transformation (1 for MKV 2 COE, else for COE 2 MKV)

% Ouput: - s, an array with the same dimensions of x, containing the transformed variables

function [s] = MKV2COE(mu, x, direction)
    
    if (direction)
        % Pre-allocation 
        s = zeros(7, size(x,2));

        % Eccentricity
        ev = x(4:6,:);                                  % Eccentricity vector
        s(2,:) = sqrt( dot(ev, ev, 1) );                % Eccentricity norm

        % Semimajor axis 
        h = x(1:3,:);                                   % Angular momentum vector
        s(1,:) = dot(h,h,1) ./ (mu * (1 - s(2,:).^2));  % Semimajor axis

        % Semilatus rectum 
        s(7,:) = s(1,:) .* (1 - s(2,:).^2);

        % Euler matrix
        k = h./ sqrt( dot(h, h, 1) );
        
        I = ev ./ s(2,:);
        K = repmat([0;0;1], 1, size(x,2));              % Inertial Z axis
        n = cross(K, s(1:3,:));                         % Node vector
        I(:, s(2,:) == 0) = n(:, s(2,:) == 0);          % Circular orbit singularity

        j = cross(k,I); 

        Q = reshape([I; j; k], 3, []).';
        Omega = atan2(Q(3:3:end,1),-Q(3:3:end,2));             % RAAN
        omega = atan2(Q(1:3:end,3),Q(2:3:end,3));              % Argument of perigee
        i = acos(Q(3:3:end,3));                                % Inclination

        s(3,:) = Omega.'; 
        s(4,:) = i.'; 
        s(5,:) = omega.'; 
        s(6,:) = x(7,:) - Omega.' - omega.';                   % Mean anomaly
        
        % Circular orbit singularity 
        idx = s(2,:) == 0;
        s(6, idx) = x(7,idx) - Omega(idx);

        % Equatorial singularity
        
    else
        % Pre-allocation 
        e = zeros(3, size(x,2));        % Eccentricity vector
        h = zeros(3, size(x,2));        % Angular momentum vector
    
        % Eccentricity vector 
        e(1,:) = x(2,:);
    
        % Angular momentum
        h(3,:) = sqrt(mu * x(1,:) .* (1 - x(2,:).^2));
    
        % Compute Euler matrix 
        for i = 1:size(x,2)
            Q = Astrodynamics.euler_matrix(x(:,i)).';
            h(:,i) = Q * h(:,i);
            e(:,i) = Q * e(:,i);
        end
    
        % Mean longitude 
        l = x(3,:) + x(5,:) + x(6,:);
    
        % Final result 
        s = [h; e; l];
    end
end