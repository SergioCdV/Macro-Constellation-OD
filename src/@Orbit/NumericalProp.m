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
    AuxOrbit = obj.ChangeStateFormat('Cartesian'); 
    [~, AuxEvolution] = ode45(@(t,s)dynamics(obj.mu, obj.mus, obj.mul, obj.J2, obj.Re, s, obj.PropagatedEpoch, t), tspan, AuxOrbit.ElementSet, obj.IntegrationOptions);
    AuxOrbit.StateEvolution = [tspan.' AuxEvolution]; 
    AuxEvolution = AuxOrbit.ChangeStateFormat(obj.ElementType).StateEvolution(:,2:end);
end

%% Auxiliary functions 
function [dS] = dynamics(mu, mu_s, mu_l, J2, Re, s, InitialEpoch, tspan)
   % State variables
   r = s(1:3);      % ECI position coordinates
   v = s(4:6);      % ECI velocity coordinates 

   % J2 Dynamics
   R = norm(r);
   gamma = -mu*r/R^3.*[1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); 1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); 1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-3)];

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
    JD = juliandate(T0+seconds(deltaT)); 

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
end