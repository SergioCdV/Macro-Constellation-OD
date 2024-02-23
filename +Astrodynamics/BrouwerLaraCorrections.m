
function [S] = BrouwerLaraCorrections(J2, J3, H, L, direction)
    % Perturbations parameters in the expansion of the Hamiltonian
    epsilon(1) = -J2/4;              % 1/4 J2 * 1       (non-dimensional Re)
    epsilon(2) = +(J3/J2) / 2;       % 1/2 (J3/J2) * 1  (non-dimensional Re)

    if (~exist('direction', 'var'))
        direction = true;
    end

    if (~direction)
        epsilon = -1 * epsilon;
    end

    % Mean to long-period transformation
    L_long = Mean2Long(epsilon(2), L);

    % Long-period to osculating transformation
    S = Long2Osc(H, epsilon(1), L_long);  
end

%% Auxiliary functions 
% Mean to long-period transformation (Brouwer)
function [Do] = Mean2Long(epsilon, L)
    psi = L(1,:);           
    chi = L(2,:);       
    xi = L(3,:);       
    r = L(4,:);           
    R = L(5,:);           
    Theta = L(6,:);             % Angular momentum       

    % Auxiliary functions 
    p = Theta.^2;               % Semilatus rectum
    k = p ./ r - 1;
    sigma = p .*R ./ Theta;
    eps = epsilon ./ p;         % 1/2 (J3/J2)(Re/p)
    c = sqrt(1-chi.^2-xi.^2);

    % Transformation
    psi_l = psi + eps * (2*xi + (k .* xi - c .* chi .* sigma) ./ (1+c));
    chi_l = chi + eps * (2*xi.^2 + k .* (1-chi.^2));
    xi_l = xi - eps * (c.^2 .* sigma + (2+k) .* chi .* xi);                                              
    r_l = r + eps * chi .* p;
    R_l = R + eps * (1+k) .* xi .* Theta ./ r;
    Theta_l = Theta + eps * (k .* chi - sigma .* xi) .* Theta;

    % Output
    Do = [psi_l; chi_l; xi_l; r_l; R_l; Theta_l];
end

% Long-period to osculating transformation (Brouwer)
function [Lo] = Long2Osc(H, epsilon, L)
    psi = L(1,:);           
    chi = L(2,:);       
    xi = L(3,:);       
    r = L(4,:);           
    R = L(5,:);           
    Theta = L(6,:);              % Angular momentum

    p = Theta.^2;                % Semilatus rectum
    c = H./Theta;                % Cosine of the inclination
    epsilon = epsilon./p.^2;     % J2/p^2

    % Transformation
    psi_o = psi - epsilon * (1+7*c)./(1+c) .* chi .* xi;
    chi_o = chi - epsilon * (xi.^2-3*c.^2) .* chi;
    xi_o = xi + epsilon * (chi.^2-3*c.^2) .* xi;
    r_o = r + 2 * epsilon * r .* (chi.^2-2+5*c.^2);                                                      
    R_o = R + 4 * epsilon * Theta./r .* chi .* xi;
    Theta_o = Theta + 3 * epsilon * Theta .* (chi.^2-xi.^2);                                                 

    % Output
    Lo = [psi_o; chi_o; xi_o; r_o; R_o; Theta_o];
end