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

    % Constants of the dynamics 
    p = dot(h,h,1);                     % Semilatus rectum
    h_norm = sqrt(p);                   % Angular momentum
    uH = h ./ h_norm;                   % Angular momentum unit vector
    e_norm = dot(e,e,1);                % Orbital eccentricity
    eta = sqrt(1-e_norm);               % Eccentricity function
    a = p ./ (1-e_norm);                % Semimajor axis 
    n = sqrt(1./a.^3);                  % Mean motion
    zeta = Keci.' * uH;                 % Cosine of the inclination

    % Jacobian 
    J = zeros(size(s,1), size(s,1)); 
    Dyad = eye(3) / h_norm - h*h.'/h_norm^3;

    dn_p = -3 * h_norm^4 * eta^3 * uH.';
    dn_e = -3 * h_norm^3 * eta * e.';
    dnp_h = +3*J2/(2*p^2) * (dn_p - 4 * n / h_norm * uH.');
    dnp_e = +3*J2/(2*p^2) * dn_e;

    J(1:3,1:3) = - 3*n*J2 / (2*p^2) * ( zeta * QuaternionAlgebra.hat_map(Keci) - QuaternionAlgebra.hat_map(h) * Dyad );
    J(1:3,1:3) = J(1:3,1:3) - zeta .* cross(h,  Keci) * dnp_h;
    J(1:3,4:6) = - zeta .* cross(uH,  Keci) .* h_norm * dnp_e;
    
    J(4:6,1:3) = + 3*n*J2 / (4*p^2) * (1-5*zeta^2 + 2) * QuaternionAlgebra.hat_map( e ) * Dyad;
    J(4:6,4:6) = - 3*n*J2 / (4*p^2) * QuaternionAlgebra.hat_map( (1-5*zeta^2)*uH + 2 * zeta * Keci );

    J(4:6,1:3) = J(4:6,1:3) - cross( (1-5*zeta.^2).*uH + 2 * zeta .* Keci, e) * dnp_h / 2;
    J(4:6,4:6) = J(4:6,4:6) - cross( (1-5*zeta.^2).*uH + 2 * zeta .* Keci, e) * dnp_e / 2;
    
    J(7,1:3) = +3*n*J2 / (4*p^2) * (6 * eta * zeta + 10 * zeta - 4) * Keci.' * Dyad;
    J(7,1:3) = J(7,1:3) + dn_p + (eta .* (3*zeta.^2-1) + 5 * zeta.^2 - 2 * zeta - 1) * dnp_h / 2;

    J(7,4:6) = -3*n*J2 / (4*p^2) * (3*zeta^2-1) * (e.' / eta);
    J(7,4:6) = J(7,4:6) + (1 + (eta .* (3*zeta.^2-1) + 5 * zeta.^2 - 2 * zeta - 1) ) * dnp_e / 2;
end