%% Constellation macro-orbit determination 
% Date: 07/02/2024
% Author: Sergio Cuevas del Valle

%% TLE to classical orbital elements %%
% This file contains the function to change from the state vector to the
% classical orbital elements bu un-kozai the mean motion

function [elements] = TLE2COE(xke, J2, s, direction)    
    % Initialization of the solution
    elements = s;

    if (direction)
        % Un-kozai the mean motion 
        n_unkozai = elements(1,:);
        a1 = (xke ./ n_unkozai) .^ (2/3);
        k2 = 0.5 * J2;
        K = 1.5 * k2 * (3.0 * cos(elements(4,:)).^2 - 1.0) ./ (1-elements(2,:).^2).^(3/2);
        d1 = K ./ a1.^2;
        a0 = a1 .* (1.0 - d1 / 3 - d1.^2 - 134/81 * d1.^3);
        del = K ./ a0.^2;
        n0 = n_unkozai ./ (1.0 + del);
        elements(1,:) = (xke ./ n0) .^ (2/3);
    else
    end
end