%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Keplerian dynamics %% 
% Function to propagate an orbit under Keplerian dynamics

% Inputs: 

% Outputs: 

function [AuxEvolution] = KeplerDynamics(obj, tspan)
    % Switch the Keplerian propagation 
    switch (obj.ElementType)
        case 'ECI'
            AuxElements = obj.ECI2COE(obj.ElementSet, true);
            AuxEvolution = COEK_dynamics(obj.mu, AuxElements, obj.PropagatedEpoch, tspan);
            AuxEvolution = obj.ECI2COE(AuxEvolution(:,1:obj.m), false);
        case 'Polar'
             AuxElements = obj.POL2COE(obj.ElementSet, true);
            AuxEvolution = COEK_dynamics(obj.mu, AuxElements, obj.PropagatedEpoch, tspan);
            AuxEvolution = obj.POL2COE(AuxEvolution(:,1:obj.m), false);  
        case 'COE'
            AuxEvolution = COEK_dynamics(obj.mu, obj.ElementSet, obj.PropagatedEpoch, tspan);
        case 'MOE'
            AuxEvolution = MOEK_dynamics(obj, obj.ElementSet, obj.PropagatedEpoch, tspan);
        case 'KS'
            AuxEvolution = KSK_dynamics(obj.mu, obj.ElementSet, obj.PropagatedEpoch, tspan);
    end
end

%% Auxiliary variables
% Kepler dynamics in the COE representation 
function [S] = COEK_dynamics(mu, s, InitialEpoch, tspan)
    % Preallocation
    S = zeros(length(tspan),7); 

    % Frozen dynamics
    S(:,[1:5 7]) = repmat(s(1, [1:5, 7]), length(tspan), 1);

    % Mean anomaly 
    S(:,6) = sqrt(mu/s(1)^3)*(tspan).';
end

% Kepler dynamics in the MOE representation 
function [S] = MOEK_dynamics(obj, s, InitialEpoch, tspan)
    % Constants 
    mu = obj.mu; 

    % Preallocation
    S = zeros(length(tspan),6); 

    % Frozen dynamics
    S(:,[1:5 6]) = repmat(s, length(tspan), 1);

    % Anomaly dynamics 
    dL = @(t,theta)(sqrt(mu*s(1))*((1+s(2)*cos(theta)+s(3)*cos(theta))/s(1))^2);
    [~, S(:,6)] = ode45(@(t,s)dL(t,s), tspan, s(6), obj.IntegrationOptions);
end

% Kepler dynamics in the KS representation 
function [S] = KSK_dynamics(mu, s, InitialEpoch, tspan)
    % Deprecated 
    error('KS propagation is not currently supported.');

    % Constants of motion 
    E = 4*(-mu/2+dot(s(5:8),s(5:8)))/dot(s(1:4),s(1:4));        % Orbital energy 
    alpha = -2/mu*E;                                            % Inverse of the semimajor axis
    W = sqrt(mu*alpha/4);                                       % KS oscillation frequency
    
    % State trajectory
    ds = 0:1e-3:1e3;
    S = [s(1:4).*cos(W*ds.')+s(5:8)/W.*sin(W*ds.') -W*s(1:4).*sin(W*ds.')+s(5:8).*cos(W*ds.')];

    t = cumtrapz(ds, dot(S(:,1:4),S(:,1:4),2));

    diff = abs(t-tspan(end));
    [~, index] = sort(diff); 
    S = S(1:index(1),:);
    t = t(1:index(1));
end

