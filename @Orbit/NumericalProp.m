%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Numerical propagator %% 
% Function to propagate an orbit under Keplerian dynamics, J2 and lunisolar
% perturbations in Cartesian coordinates

% Inputs: 

% Outputs:

function [AuxEvolution] = NumericalProp(obj, tspan)
    % Change to Cartesian coordinates 
    AuxOrbit = obj.ChangeStateFormat('ECI').Normalize(false, 1); 
    if (obj.Normalized)
        tspan = tspan * obj.Tc;
    end

    % Integration
    [~, AuxEvolution] = ode45(@(t,s)dynamics(AuxOrbit.mu, AuxOrbit.mus, AuxOrbit.mul, AuxOrbit.J2, AuxOrbit.Re, s, AuxOrbit.PropagatedEpoch, t), tspan, AuxOrbit.ElementSet, AuxOrbit.IntegrationOptions);

    % Formatting
    AuxOrbit.StateEvolution = [tspan.' AuxEvolution]; 
    AuxEvolution = AuxOrbit.ChangeStateFormat(obj.ElementType);
    if (obj.Normalized)
        AuxEvolution = AuxEvolution.Normalize(true, obj.Lc);
    end
    AuxEvolution = AuxEvolution.StateEvolution(:,2:end);
end

%% Auxiliary functions 
function [dS] = dynamics(mu, mu_s, mu_l, J2, Re, s, InitialEpoch, tspan)
   % State variables
   r = s(1:3);      % ECI position coordinates
   v = s(4:6);      % ECI velocity coordinates 

   % J2 Dynamics
   R = norm(r);
   gamma = -mu*r/R^3.*[1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); ...
                       1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); ...
                       1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-3)];

   % Lunisolar perturbations
   s = sun_ephemerides(InitialEpoch, tspan);
   l = moon_ephemerides(InitialEpoch, tspan);
   gamma_ls = mu_s*((s-r)/norm(s-r)^3-s/norm(s)^3)+mu_l*((l-r)/norm(l-r)^3-l/norm(l)^3);

   % Total dynamics
   dS = [v; gamma+gamma_ls];
end

% Sun ephemeris
function [s] = sun_ephemerides(T0, deltaT)
    % Current time 
    JD = T0+deltaT/86400; 

    % Orbital elements
    T = (JD-2451545.0)/36525.0;
    M = deg2rad(357.5277233)+deg2rad(35999.05034)*T;
    Omega = deg2rad(280.460)+deg2rad(36000.771)*T;

    % Geocentric coordinates 
    eps = deg2rad(23.43929111)-deg2rad(0.0130042)*T;
    lambda = Omega+deg2rad(6892/3600)*sin(M)+deg2rad(72/3600)*sin(2*M);
    r = (149.619-2.499*cos(M)-0.021*cos(2*M))*1e9;
   
    s = r*[cos(lambda); sin(lambda)*cos(eps); sin(lambda)*sin(eps)];
end

% Moon ephemeris
function [l] = moon_ephemerides(T0, deltaT)
    % Current time 
    JD = T0+deltaT/86400; 

    % Orbital elements
    T = (JD-2451545.0)/36525.0;
    L0 = deg2rad( 218.31617 + 481267.88088 * T -1.3972 * T );
    l = deg2rad( 134.96292 + 477198.86753 * T );
    lp = deg2rad( 357.52543 + 35999.04944 * T );
    F = deg2rad( 93.27283 + 483202.01873 * T );
    D = deg2rad( 297.85027 + 445267.11135 * T );

    % Longitude of the Moon 
    lambda = L0+(22640*sin(l)+769*sin(2*l)-4586*sin(l-2*D)+2370*sin(2*D)-668*sin(lp)-412*sin(2*F)-212*sin(2*l-2*D)-206*sin(l+lp-2*D)+192*sin(l+2*D)-165*sin(lp-2*D)+148*sin(l-lp)-125*sin(D)-110*sin(l+lp)-55*sin(2*F-2*D)) / 3600;

    % Latitude of the Moon 
    beta = (18520*sin(F+lambda-L0+412/3600*sin(2*F)+541/3600*sin(lp))-526*sin(F-2*D)+44*sin(l+F-2*D)-31*sin(F-2*D-l)-25*sin(F-2*l)-23*sin(F-2*D+lp)+21*sin(F-l)+11*sin(F-2*D-lp))/3600;

    % Distance to the Moon 
    Rm = 385000-20905*cos(l)-3699*cos(2*D-l)-2956*cos(2*D)-570*cos(2*l)+246*cos(2*l-2*D)-205*cos(lp-2*D)-171*cos(l+2*D)-152*cos(l+lp-2*D); 
    Rm = Rm * 1e3;

    % Final position vector in Geocentric coordinates
    eps = deg2rad(23.43929111)-deg2rad(0.0130042)*T;
    R = [1 0 0; 0 cos(eps) -sin(eps); 0 sin(eps) cos(eps)];
    r = Rm*[cos(lambda)*cos(beta); sin(lambda)*cos(beta); sin(beta)];
    l = R*r;
end