
function [S] = Lara_solution(epsilon, D)
    % Delaunay to Lara variables 
    epsilon = epsilon/4;
    L = Astrodynamics.Delaunay2Lara(D, true);

    % Mean to long-period transformation
    Ll = Mean2Long(-epsilon, L);

    % Long-period to osculating transformation
    Lo = Long2Osc(epsilon, Ll);

    % Lara variables to ECI 
    S = Astrodynamics.Lara2ECI(Lo, true);
end

%% Auxiliary functions 
% Mean to long-period transformation (Brouwer)
function [Do] = Mean2Long(epsilon, L)
    psi = L(1);           
    chi = L(2);       
    xi = L(3);       
    r = L(4);           
    R = L(5);           
    Theta = L(6);         

    % Auxiliary functions 
    p = Theta^2;
    k = -1+p/r;
    sigma = p*R/Theta;
    epsilon = epsilon/p^2;

    % Transformation
    psi_l = psi;                                                     % no J3
    chi_l = chi + 7/8 * epsilon * (xi*(k^2-sigma^2)+2*k*sigma*chi);
    xi_l = xi -7/8 * epsilon * (chi*(k^2-sigma^2)-2*k*sigma*xi);
    r_l = r;                                                         % no J3
    R_l = R;                                                         % no J3
    Theta_l = Theta;                                                 % no J3

    % Output
    Do = [psi_l, chi_l, xi_l, r_l, R_l, Theta_l].';
%     Do = [psi, chi, xi, r, R, Theta].';
end

% Long-period to osculating transformation (Brouwer)
function [Lo] = Long2Osc(epsilon, L)
    psi = L(1);           
    chi = L(2);       
    xi = L(3);       
    r = L(4);           
    R = L(5);           
    Theta = L(6);     

    p = Theta^2;
    c = sqrt(1-xi^2-chi^2);
    epsilon = epsilon/p^2; 

    % Transformation
    psi_o = psi - epsilon * (1+7*c)/(1+c) * chi * xi;
    chi_o = chi - epsilon * (xi^2-3*c^2) * chi;
    xi_o = xi + epsilon * (chi^2-3*c^2) * xi;
    r_o = r + 2 * epsilon * r * (chi^2-2+5*c^2);                                                      
    R_o = R + 4 * epsilon * Theta/r * chi * xi;
    Theta_o = Theta + 3 * epsilon * Theta * (chi^2-xi*2);                                                 

    % Output
    Lo = [psi_o, chi_o, xi_o, r_o, R_o, Theta_o].';
end