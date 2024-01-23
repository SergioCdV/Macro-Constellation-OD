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
    % Reshape 
    s = reshape(s, 7, []); 

    % State vector partition 
    h = s(1:3,:);                       % Angular momentum vector 
    e = s(4:6,:);                       % Eccentricity vector 
    l = s(7,:);                         % Longitude

    % Constants of the dynamics 
    p = dot(h,h,1);                     % Semilatus rectum
    h_norm = sqrt(p);                   % Angular momentum
    uH = h ./ h_norm;                   % Angular momentum unit vector
    e_norm = sqrt(dot(e,e,1));          % Orbital eccentricity
    a = p.^2 ./ (1-e_norm.^2);          % Semimajor axis 
    eta = sqrt(1-e_norm.^2);            % Eccentricity function
    n = a.^(-3/2);                      % Mean motion
    zeta = Keci.' * uH;                 % Cosine of the inclination

    % Dynamics
    ds(1:3,:) = + 3*n*J2 ./ (4*p.^2) .* ( QuaternionAlgebra.hat_map(Keci) * (zeta.*uH) );
    ds(4:6,:) = - 3*n*J2 ./ (4*p.^2) .* cross( (1-5*zeta.^2).*uH + 2 * zeta .* repmat(Keci, 1, length(zeta)), e);
    ds(7,:) = n + 3*n*J2 ./ (4*p.^2) .* (eta .* (3*zeta.^2-1) + 5 * zeta.^2 - 2 * zeta - 1);

    ds = reshape(ds, [], 1);
end