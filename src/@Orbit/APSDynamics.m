%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Artificial Satellite Problem dynamics %% 
% Function to propagate an orbit under Keplerian dynamics and J2 perturbations

% Inputs: 

% Outputs: 

function [AuxEvolution] = APSDynamics(obj, tspan, model)
    % Branch the dynamics
    switch (model)
        case 0
            % Switch the Keplerian propagation
            AuxOrbit = obj.ChangeStateFormat('Cartesian');
        
            % Integration of the osculating problem
            [~, AuxEvolution] = ode45(@(t,s)APSO_dynamics(obj.mu, obj.J2, obj.Re, s, obj.PropagatedEpoch, t), tspan, AuxOrbit.ElementSet, obj.IntegrationOptions);
            AuxOrbit.StateEvolution = [tspan.' AuxEvolution]; 
            AuxEvolution = AuxOrbit.ChangeStateFormat(obj.ElementType).StateEvolution(:,2:end);

        case 1
            % Switch the Keplerian propagation
            AuxElements = obj.ChangeStateFormat('COE').ElementSet;
        
            % Integration of the mean problem
            AuxEvolution = APSM_dynamics(obj.mu, obj.J2, obj.Re, AuxElements, obj.PropagatedEpoch, tspan);

        case 2
            error('SGP4 is not currently supported.');
        otherwise
    end
end

%% Auxiliary variables
% J2 osculating dynamics 
function [dS] = APSO_dynamics(mu, J2, Re, s, InitialEpoch, tspan)
   r = s(1:3);      % ECI position coordinates
   v = s(4:6);      % ECI velocity coordinates 

   % Dynamics
   R = norm(r);
   gamma = -mu*r/R^3.*[1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); 1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); 1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-3)];
   dS = [v; gamma];
end

% J2 mean dynamics 
function [AuxEvolution] = APSM_dynamics(mu, J2, Re, s, InitialEpoch, tspan)
    % Constants 
    a = s(1);            % Semimajor axis 
    e = s(2);            % Eccentricity 
    i = s(4);            % Inclination
    n = sqrt(mu/a^3);    % Mean motion
    p = a*(1-e^2);       % Semilatus rectum

    % Kozai's method
    dS = [-3/2*n*J2*(Re/p)^2*cos(i); ...
           3/2*n*J2*(Re/p)^2*(2-5/2*sin(i)^2); ...
          n*(1+3/2*J2*(Re/p)^2*sqrt(1-e^2)*(1-3/2*sin(i)^2))];      
    
    % Integration 
    AuxEvolution(:, [1 2 4 7]) = repmat([a e i s(7)], length(tspan), 1);
    AuxEvolution(:, [3 5 6]) = repmat(s([3 5 6]), length(tspan), 1)+(dS.*tspan).';
end

% SGP4 
function [S] = SGP4_dynamics(mu, s, InitialEpoch, tspan)
    S = s;
end

