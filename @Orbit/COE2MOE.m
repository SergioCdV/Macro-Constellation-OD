%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Classical orbital elements to equinoctial %% 
% Function to transform classical orbital elements into Cartesian or
% viceversa

% Inputs: - vector s, the orbital state vector to be transformed 
%         - boolean direction, 0 for the COE to equinoctional
%           transformation and 1 for the viceversa case

% Outputs: - vector S, the transformed orbital state vector

function [S] = COE2MOE(obj, s, direction)
    % Constants
    mu = obj.mu; 
    
    % Sanity check on the s dimensions 
    if (size(s,1) ~= 6)
        lastwarn('Input orbital element set is not valid')
    end

    % Switch directions 
    if (direction)
        % Preallocation 
        S = zeros(size(s,1),6);

        % Transformation
        for i = 1:size(s,1)
            S(i,1) = s(i,1)*(1-s(i,2)^2);            
            S(i,2) = s(i,2)*cos(s(i,5)+s(i,3));
            S(i,3) = s(i,2)*sin(s(i,5)+s(i,3));
            S(i,4) = tan(s(i,4)/2)*cos(s(i,3));
            S(i,5) = tan(s(i,4)/2)*sin(s(i,3));
            S(i,6) = s(i,6)+s(i,5)+s(i,3);
        end
    else
        % Preallocation 
        S = zeros(size(s,1),7);

        for i = 1:size(s,1)
            S(i,1) = s(i,1)/(1-s(i,2)^2-s(i,3)^2);                                    % Semimajor axis
            S(i,2) = sqrt(s(i,2)^2+s(i,3)^2);                                         % Eccentricity
            S(i,3) = atan2(s(i,5),s(i,4));                                            % RAAN
            S(i,4) = atan2(2*sqrt(s(i,4)^2+s(i,5)^2), 1-s(i,4)^2-s(i,5)^2);           % Inclination
            S(i,5) = atan2(s(i,3)*s(i,4)-s(i,2)*s(i,5),s(i,2)*s(i,4)+s(i,3)*s(i,5));  % Argument of periapsis
            S(i,6) = s(i,6)-(S(i,5)+S(i,3));                                          % True anomaly
            S(i,7) = s(i,1)*(1-s(i,2)^2);                                             % Semilatus rectum
        end
    end  
end