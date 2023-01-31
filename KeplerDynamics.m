%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Classical orbital elements to equinoctial %% 
% Function to transform classical orbital elements into Cartesian or
% viceversa

% Inputs: - vector s, the orbital state vector to be transformed 
%         - boolean direction, 0 for the COE to equinoctional
%           transformation and 1 for the viceversa case

% Outputs: - vector S, the transformed orbital state vector

function [obj] = KeplerDynamics(obj)
    % Switch the Keplerian propagation 
    switch (obj.ElementType)
        case 'Cartesian'
            AuxElements = ECI2COE(obj.mu, obj.ElementSet, true);
            AuxEvolution = COEK_dynamics(AuxElements, obj.PropagatedEpoch, obj.CurrentEpoch);
        case 'COE'
            AuxEvolution = COEK_dynamics(obj.ElementSet, obj.PropagatedEpoch, obj.CurrentEpoch);
        case 'MOE'
            AuxEvolution = MOEK_dynamics(obj.ElementSet, obj.PropagatedEpoch, obj.CurrentEpoch);
        case 'KS'
            AuxEvolution = KSK_dynamics(obj.ElementSet, obj.PropagatedEpoch, obj.CurrentEpoch);
    end

    % Tspan 
    tspan = obj.PropagatedEpoch:obj.TimeStep:obj.CurrentEpoch;

    % Update states
    obj.StateEvolution = [obj.StateEvolution; [tspan.' AuxEvolution]];
    obj.ElementSet = AuxEvolution(end,:);
end

%% Auxiliary variables
% Kepler dynamics in the COE representation 
function [S] = COEK_dynamics(s, PastEpoch, CurrentEpoch)
end

% Kepler dynamics in the MOE representation 
function [S] = MOEK_dynamics(s, PastEpoch, CurrentEpoch)
end

% Kepler dynamics in the KS representation 
function [S] = KSK_dynamics(s, PastEpoch, CurrentEpoch)
end