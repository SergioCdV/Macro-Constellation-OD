
function [Do] = Brouwer_solution(epsilon, D)
    % Mean to long-period transformation
    Dl = Mean2Long(epsilon, D);

    % Long-period to osculating transformation
    Do = Long2Osc(epsilon, Dl);
end

%% Auxiliary functions 
% Mean to long-period transformation (Brouwer)
function [Do] = Mean2Long(epsilon, D)
    M = D(1);           % Mean anomaly
    Omega = D(2);       % RAAN
    omega = D(3);       % AoP
    L = D(4);           % Delaunay action
    G = D(5);           % Angular momentum
    H = D(6);           % Polar angular momentum

    % Transformation
    Ll = L; 
    Gl = G + epsilon/(16*G^3) * (1 - G^2/L^2) * (1-16*H^2/G^2+15*H^4/G^4) * cos(2*omega)/(1-5*H^2/G^2);
    Hl = H;
    Ml = M - epsilon/(16*G^4) * (G/L)^3 * (1-16*H^2/G^2+15*H^4/G^4) * sin(2*omega)/(1-5*H^2/G^2);
    omegal = omega + epsilon/(32*G^4)*1/(1-5*H^2/G^2)*((3-G^2/L^2)*(1-16*H^2/G^2+15*H^4/G^4)-2*(H/G)^2*(1-G^2/L^2)*(11+25*H^2/G^2+200*(H^4/G^4)/(1-5*H^2/G^2)))*sin(2*omega);
    Omegal = Omega + epsilon/(16*G^4) * (H/G)*(1-(G/L)^2)/(1-5*(H/G)^2)*(11+25*H^2/G^2+200*(H^4/G^4)/(1-5*H^2/G^2)*sin(2*Omega));

    % Output
    Do = [Ml Omegal omegal Ll Gl Hl].';
end

% Long-period to osculating transformation (Brouwer)
function [Do] = Long2Osc(epsilon, D)
    M = D(1);           % Mean anomaly
    Omega = D(2);       % RAAN
    omega = D(3);       % AoP
    L = D(4);           % Delaunay action
    G = D(5);           % Angular momentum
    H = D(6);           % Polar angular momentum

    % Keplerian functions
    a = L^2;
    e = real(sqrt(1-(G/L)^2));
    cos_i = (1-(H/G)^2); 

    % Solve for the true anomaly and the radial magnitude
    nu = Astrodynamics.KeplerSolver(e,M);
    R = G./(1+e*cos(nu));
  
    % Transformation
    Lo = L - epsilon/(4*L^3) * ((3*(H/G)^2-1) * ((a/R)^3-(L/G)^3) + 3 * cos_i * (a/R)^3 * cos(2*nu+2*omega)); 
    Go = G - 3*epsilon/(4*G^3) * cos_i * (cos(2*nu+2*omega) + e*cos(nu+2*omega) + e/3*cos(3*nu+2*omega));
    Ho = H;
    Mo = M + epsilon/(8*e*L^4) * (L/G) * (2*(3*(H/G)^2-1)*(1+a/R+(a*G)^2/(R*L)^2)*sin(nu) + 3*cos_i*((1-a/R-(a*G)^2/(R*L)^2)*sin(nu+2*omega) + (1/3+a/R+(a*G)^2/(R*L)^2)*sin(3*nu+2*omega)));
    omegao = omega - epsilon/(8*e*L^4) * (L/G)^2 * (2*(3*(H/G)^2-1)*(1+a/R+(a*G)^2/(R*L)^2)*sin(nu) + 3*cos_i*((1-a/R-(a*G)^2/(R*L)^2)*sin(nu+2*omega) + (1/3+a/R+(a*G)^2/(R*L)^2)*sin(3*nu+2*omega))) ...
             - 3*epsilon/(8*G^4)*(2*(5*(H/G)^2-1))*(nu-M+e*sin(nu) + 3*(3-5*(H/G)^2)*(sin(2*nu+2*omega) + e*sin(nu+2*omega) + e/3*sin(3*nu+2*omega)));
    Omegao = Omega + 3*epsilon/(4*G^4) * (H/G) * (2*(nu-M+e*sin(nu)) -sin(2*nu+2*omega) - e*sin(nu+2*omega) - e/3*sin(3*nu+2*omega));

    % Output
    Do = [Mo Omegao omegao Lo Go Ho].';
end