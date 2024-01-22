%% Constellation macro-orbit determination %%
% Date: 22/01/2024
% Author: Sergio Cuevas del Valle

%% Dynamics of the Milankovitch's elements %% 
% Function to compute the Milankovitch's set of classical elements dynamics
% vector field

% Inputs: - scalar J2, the J2 parameter of the central body under
%           consideration
%         - vector Keci, denoting the third vector in the ECI triad of unit
%           vectors
%         - scalar t, the current epoch 
%         - vector s, the state vector to be integrated

% Outputs: - vector ds, the dynamics of the orbital set

function [ds] = milankovitch_dynamics(J2, Keci, t, s)
    % State vector partition 
    h = s(1:3,1);                       % Angular momentum vector 
    e = s(4:6,1);                       % Eccentricity vector 
    l = s(7,1);                         % Longitude

    % Constants of the dynamics 
    h_norm = norm(h);                   % Angular momentum
    uH = h / h_norm;                    % Angular momentum unit vector
    e_norm = norm(e);                   % Orbital eccentricity
    p = h_norm^2;                       % Semilatus rectum
    a = p^2 / (1-e_norm^2);             % Semimajor axis 
    eta = sqrt(1-e_norm^2);             % Eccentricity function
    n = a^(-3/2);                       % Mean motion
    zeta = dot(Keci, uH);

    % Dynamics
    ds(1:3,1) = - 3*n*J2 / (4*p^2) * cross(zeta*uH, Keci);
    ds(4:6,1) = - 3*n*J2 / (4*p^2) * cross((1-5*zeta^2)*uH + 2 * zeta * Keci, e);
    ds(7,1) = n + 3*n*J2 / (4*p^2) * (eta * (3*zeta^2-1) + 3*zeta^2 - 1);
end