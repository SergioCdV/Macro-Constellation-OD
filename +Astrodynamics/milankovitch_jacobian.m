%% Constellation macro-orbit determination %%
% Date: 22/01/2024
% Author: Sergio Cuevas del Valle

%% Jacobian of the Milankovitch's elements dyanmics %% 
% Function to compute the Jacobian of the Milankovitch's set of classical elements dynamics


% Inputs: - scalar J2, the J2 parameter of the central body under
%           consideration
%         - vector Keci, denoting the third vector in the ECI triad of unit
%           vectors
%         - vector s, the state vector to be integrated

% Outputs: - vector ds, the dynamics of the orbital set

function [J] = milankovitch_jacobian(J2, Keci, s)
    % State vector partition 
    h = s(1:3,1);                       % Angular momentum vector 
    e = s(4:6,1);                       % Eccentricity vector 
    l = s(7,1);                         % Longitude

    % Constants of the dynamics 
    h_norm = norm(h);                   % Angular momentum
    uH = h / h_norm;                    % Angular momentum unit vector
    e_norm = norm(e);                   % Orbital eccentricity
    i = e / e_norm;                     % Orbital eccentricity unit vector
    p = h_norm^2;                       % Semilatus rectum
    a = p^2 / (1-e_norm^2);             % Semimajor axis 
    eta = sqrt(1-e_norm^2);             % Eccentricity function
    n = a^(-3/2);                       % Mean motion
    zeta = dot(Keci, uH);

    % Jacobian 
    J = zeros(size(s,1), size(s,1)); 
    J(1:3,1:3) = - 3*n*J2 / (4*p^2) * ( zeta * QuaternionAlgebra.hat_map(Keci) - QuaternionAlgebra.hat_map(uH) );
    J(4:6,1:3) = + 3*n*J2 / (4*p^2) * (1-5*zeta^2 + 2) * QuaternionAlgebra.hat_map( i );
    J(4:6,4:6) = - 3*n*J2 / (4*p^2) * QuaternionAlgebra.hat_map( (1-5*zeta^2)*uH + 2 * zeta * Keci );
    J(7,1:3) = 3*n*J2 / (4*p^2) * (6 * eta * zeta + 6 * zeta);
end