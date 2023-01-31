%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% KS matrix %%
% Function to compute the operator to transform from the position space to
% the u-space

% Inputs: - vector u, the state variable

% Outputs: - array L, the KS operator

function [L] = KS_matrix(u)
    % Compute the L operator 
    L = [u(1) -u(2) -u(3) u(4); u(2) u(1) -u(4) -u(3); u(3) u(4) u(1) u(2); u(4) -u(3) u(2) -u(1)];
end