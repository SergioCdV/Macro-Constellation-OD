%% Constellation tracking
% Author: Sergio Cuevas
% Date: 13/03/2024

%% Lara's osculating to mean transformations
% The following function provides the mean-to-osculating and osculating-to-mean transformatiosn for the zonal J2/J3 Brouwer's problem
% The transformations are of first order in the eccentricity, just as Simplified General Perturbations 4 (SGP4), and make use of Lara's
% nonsingular set of variables (a torsion transformations of the Delaunay orbital elements

% Inputs: - mu, scalar, the gravitational constant of the central body
%         - Re, scalar, the mean radius of the central body
%         - J2, scalar, the unnormalized J2 zonal coefficient of the central body
%         - J3, scalar, the unnormalized J3 zonal coefficient of the central body
%         - H, scalar, the projection of the orbital motion angular momentum onto the symmetry axis of the elliosoid of reference (usually z). H is an integral of the zonal problem
%         - L, a 6 x n vector of Lara's variables, in the order (psi, chi, phi, R, r, Theta) 
%         - direction, a boolean to indicate the direction of the transformation: 1 for mean to osculating and -1 from osculating to mean

% Ouput: - S, an array with the same dimensions of L, containing the transformed Lara's variables

function [S] = BrouwerLaraCorrections(mu, Re, J2, J3, H, L, direction)
    % Perturbations parameters in the expansion of the Hamiltonian
    epsilon(1) = -J2  * Re^2 / 4;       % - 1/4 * J2 * Re^2   (direct change coefficient)
    epsilon(2) = +J3/J2 * Re / 2;       % + 1/2 (J3/J2) * Re  (direct change coefficient)

    if (~exist('direction', 'var'))
        direction = true;
    end

    if (~direction)
        epsilon = -1 * epsilon;
        
        % Osculating to long-period transformation (short period transformation)
        L_long = Long2Osc(mu, epsilon(1), H, L); 

        % Long-period to mean transformation
        S = Mean2Long(mu, epsilon(1), epsilon(2), H, L_long);

    else
        % Mean to long-period transformation
        L_long = Mean2Long(mu, epsilon(1), epsilon(2), H, L);
    
        % Long-period to osculating transformation (short period transformation)
        S = Long2Osc(mu, epsilon(1), H, L_long);  
    end

end

%% Auxiliary functions 
% Mean to long-period transformation (Brouwer)
function [Do] = Mean2Long(mu, epsilon_2, epsilon_3, H, L)
    psi = L(1,:);           
    chi = L(2,:);       
    xi = L(3,:);       
    r = L(4,:);                     % Module of the position vector
    R = L(5,:);                     % Radial velocity
    Theta = L(6,:);                 % Angular momentum       

    % Auxiliary functions 
    p = Theta.^2 / mu;              % Semilatus rectum

    eps_3 = epsilon_3 ./ p;         % 1/2 (J3/J2)(Re/p)
    eps_2 = epsilon_2 ./ p.^2;      % -1/4 (J2)(Re^2/p^2)
    
    k = p ./ r - 1;
    sigma = p .*R ./ Theta;
    eta = sqrt(1-chi.^2-xi.^2);     % Eccentricity function

    c = H ./ Theta;                 % Cosine of the inclination
    s = sqrt(1 - c.^2);             % Sine of the inclination

    % Polynomials of the inclination
    q0 = (1 - 15 * c.^2) .* (1 - 5 * c.^2);
    q2 = s.^2 .* q0;
    q5 = c.^2 .* (11 - 30 * c.^2 + 75 * c.^4);
    q6 = q5 ./ c;
    q7 = 0.25 * (1 + 3 * c.^2 - 5 * c.^4 + 225 * c.^6);
    q8 = 0.25 * (1 - 45 * c.^2 + 195 * c.^4 - 375 * c.^6);
    q9 = 0.25 * (1 + 75 * c.^4);
    q10 = 0.25 * (1 - 40 * c.^2 + 75 * c.^4);
    q11 = 2 * c.^2 .* (6 - 25 * c.^2 + 75 * c.^4);
    q12 = 10 * c.^2;
    q13 = q0 .* (1 + c);
    q14 = 0.25 * (1 - c) .* (1 - 20 * c - 40 * c.^2 + 75 * c.^4);
    q15 = 0.25 * (1 + 23 * c - 20 * c.^2 - 80 * c.^3 + 75 * c.^4 + 225 * c.^5);

    P1 = q2 .* k + q7 .* k.^2 + q8 .* sigma.^2;
    P2 = q0 .* k + q9 .* k.^2 + q10 .* sigma.^2;
    P3 = q2 + q11 .* k;
    P4 = q0 + q12 .* k;

    % Avoid resonances for the critical inclination
    idx = (1 - 5 * c.^2) == 0;
    eps_2(1,idx) = zeros(1, sum(idx)); 
    c(1,idx) = cosd(63.39999999) * ones(1,sum(idx));

    % Transformations 
    % First angle
    dpsi = (eps_2 ./ (2 * (1 - 5 * c.^2)) .* (2 * chi .* xi .* (q13 .* k + q14 .* k.^2 + q15 .* sigma.^2) - sigma .* (xi.^2 - chi.^2) .* (q13 - q6 .* k)) ...
            + eps_3 .* ((2 + 2 * c + k) .* xi - c .* sigma .* chi) ) ./ (1 + c);
    psi_l = psi + dpsi;

%     psi_l = psi + eps_3 .* (2*xi + (k .* xi - c .* chi .* sigma) ./ (1+c));

    % Second angle
    dchi = eps_2 ./ (4 * (1 - 5 * c.^2)) .* (P1 .* chi + P2 .* (3 * xi.^2 - chi.^2) .* chi - P3 .* sigma .* xi - P4 .* sigma .* (xi.^2 - 3 * chi.^2) .* chi) ...
           + 0.5 * eps_3 .* (2 * s.^2 + (1 + c.^2) .* k + (2+k) .* (xi.^2 - chi.^2));
    chi_l = chi + dchi;

%     chi_l = chi + eps_3 .* (2*xi.^2 + k .* (1-chi.^2));

    % Third angle
    dxi = - eps_2 ./ (4 * (1 - 5 * c.^2)) .* (P1 .* xi + P2 .* (3 * chi.^2 - xi.^2) .* xi + P3 .* sigma .* chi + P4 .* sigma .* (chi.^2 - 3 * xi.^2) .* chi) - eps_3 .* (c.^2 .* sigma + (2+k) .* chi .* xi);
    xi_l = xi + dxi;

%     xi_l = xi - eps_3 .* (c.^2 .* sigma + (2+k) .* chi .* xi);      

    % Correction to the position vector
    dr = eps_2 .* (1 - 15 * c.^2) ./ (4 * (1 - 5 * c.^2)) .* (2 * sigma .* chi .* xi - k.* (chi.^2 - xi.^2)) + eps_3 .* chi;
    r_l = r + dr;

%     r_l = r + eps_3 .* chi .* p;

    % Correction to the radial velocity
    dR = Theta ./ p .* (1 + k).^2 .* (eps_3 .* xi - eps_2 .* (1 - 15 * c.^2) ./ (4 * (1 - 5 * c.^2)) .* (2 * k .* chi .* xi + sigma .* (chi.^2 - xi.^2)));
    R_l = R + dR;

%     R_l = R + eps_3 .* (1+k) .* xi .* Theta ./ r;

    % Correction to the angular momentunm
    dTheta = Theta.* (eps_2 .* (1 - 15 * c.^2) ./ (4 * (1 - 5 * c.^2)) .* (( k.^2 - sigma.^2) .* (xi.^2 - chi.^2) + 4 * k .* sigma .* xi .* chi) + eps_3 .* (k .* chi - sigma .* xi));
    Theta_l = Theta + dTheta;

%     Theta_l = Theta + eps_3 .* (k .* chi - sigma .* xi) .* Theta;

    % Output
    Do = [psi_l; chi_l; xi_l; r_l; R_l; Theta_l];
end

% Long-period to osculating transformation (Brouwer)
function [Lo] = Long2Osc(mu, eps, H, L)
    psi = L(1,:);           
    chi = L(2,:);       
    xi = L(3,:);       
    r = L(4,:);                             % Module of the position vector
    R = L(5,:);                             % Radial velocity
    Theta = L(6,:);                         % Angular momentum

    p = Theta.^2 / mu;                      % Semilatus rectum
    c = H./Theta;                           % Cosine of the inclination
    s = sqrt(1 - c.^2);                     % Sine of the inclination
    eps = eps ./ p.^2;                      % J2/p^2
    k = p ./ r - 1;             
    sigma = p .* R ./ Theta;
    e_squared = k.^2 + sigma.^2;            % Square of the eccentricity
    e = sqrt(e_squared);                    % Eccentricity
    eta = sqrt(1 - e_squared);              % Eccentricity function

    % Equation of the center 
    f = atan2(sigma, k);                                    % True anomaly
    sin_E = eta .* sin(f) ./ (1+e.*cos(f));                 % Sine of the eccentric anomaly
    cos_E = (cos(f)+e) ./ (1+e.*cos(f));                    % Cosine of the eccentric anomaly
    E = atan2(sin_E, cos_E);                                % Eccentric anoamly
    l = E - e .* sin(E);                                    % Mean anomaly
    theta = f - l;

    % Transformation
    dpsi = + eps .* ((3 + 6 * c - 15 * c.^2) .* theta + sigma .* (2 + 6 * c - 12 .* c.^2 + (1 - 3 * c.^2) .* (2 + k) ./ (1 + eta) + (2 + 4 * c) ./ (1 + c) .* (xi.^2 - chi.^2)) ...
                        - (1 + 7 * c + 4 * (1 + 3 * c) .* k) ./ (1 + c) .* chi .* xi);
    psi_o = psi + dpsi;

%     psi_o = psi - eps .* (1+7*c)./(1+c) .* chi .* xi;

    % Second angle 
    dchi = +eps .* (sigma .* (4 * xi.^2 - 12 * c.^2 + (1 - 3 * c.^2) .* (2 + k) ./ (1 + eta) ) .* xi - ((1 + 4 * k) .* xi.^2 - (3 + 4 * k) .* c.^2) .* chi + 3 * (1 - 5 * c.^2) .* theta .* xi);
    chi_o = chi + dchi;

%     chi_o = chi - eps .* (xi.^2-3*c.^2) .* chi; 

    % Third angle
    dxi = -eps .* (sigma .* (4 * xi.^2 - 8 * c.^2 + (1 - 3 * c.^2) .* (2 + k) ./ (1 + eta) ) .* chi - ((1 + 4 * k) .* chi.^2 - (3 + 4 * k) .* c.^2) .* xi + 3 * (1 - 5 * c.^2) .* theta .* chi);
    xi_o = xi + dxi;

%     xi_o = xi + eps .* (chi.^2-3*c.^2) .* xi;  

    % Corrections to the radial position
    dr = eps .* p .* (chi.^2 - xi.^2 + (1 + k ./ (1+eta) + 2*eta ./ (1+k)) .* (2 - 3*s.^2));
    r_o = r + dr;

    % r_o = r + 2 * eps .* r .* (chi.^2-2+5*c.^2); 

    % Corrections to the radial velocity
    dR = eps .* Theta./p .* (4 * (1+k).^2 .* chi .* xi - sigma .* (eta + (1+k).^2 ./ (1 + eta) ) .* (2 - 3 * s.^2) );
    R_o = R + dR;

%     R_o = R + 4 * eps .* Theta./r .* chi .* xi;

    % Correction to the angular momentum (complete)
    dTheta = eps .* Theta .* ( (3 + 4 * k) .* (chi.^2 - xi.^2) - 4 * sigma .* chi .* xi );
    Theta_o = Theta + dTheta; 

%     Theta_o = Theta + 3 * eps .* Theta .* (chi.^2-xi.^2);                                                

    % Output
    Lo = [psi_o; chi_o; xi_o; r_o; R_o; Theta_o];
end