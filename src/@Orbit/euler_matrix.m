%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Euler matrix %%
% This file drives a the function to compute the classical Euler rotation matrix (ZXZ sequence) from the orbital elements (RAAN, i, omega)

% Inputs: - vector elements, the mean classical Euler elements

% Outputs: - matrix Q, the rotation matrix from the inertial to the perifocal frame

function [Q] = euler_matrix(elements)
    % Elements of interest 
    RAAN = elements(3); 
    i = elements(4); 
    omega = elements(5); 
    
    % Compute the rotation matrix (Euler sequence ZXZ)
    Q1 = [cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];
    Q2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    Q3 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
    Q = Q3*Q2*Q1;
end