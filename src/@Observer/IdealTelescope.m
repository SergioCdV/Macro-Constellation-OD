%% Constellation macro-orbit determination 
% Date: 06/02/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of an IdealTelescope object

classdef IdealTelescope
    properties
        Type;                   % Class of observer

        State;                  % State of the observer
        Measurements = {};      % Measurements acquired by the observer

        InstrumentParams;       % Instrument parameters

        InitialEpoch;           % Initial epoch
        PropagatedEpoch;        % Propgated epoch
        CurrentEpoch;           % Current epoch

        Likelihood;             % Instrument likelihood function
        PD;                     % Probability of detection

        Dynamics;               % Dynamics of the observer
        ObservationModel;       % Observation function

        Tspan;                  % Considered time span
    end

    % Private properties
    properties (Access = private)
        StateDim = 3;           % Dimension of the state vector
        MeasDim = 3;            % Dimension of the measurement vector
    end

    methods 
        % Constructor
        function [obj] = IdealTelescope()
        end

        % Add probability of detection 
        function [obj] = probability_detection(myPD)
            obj.PD = myPD;
        end

        % ObservationModel 
        function [timestamp, meas, StateEvolution] = ObservationModel(Orbit, Tspan)
            % Propagate the observer 
            StateEvolution = obj.Dynamics(obj.State, Tspan); 

            % Propagate the orbit 
            Orbit = Orbit.SetCurrentEpoch(Tspan(end));
            Orbit = Orbit.Propagate();
            Orbit = Orbit.ChangeStateFormat('Cartesian');

            % Restrict the orbit state evolution
            diff = abs(Tspan(1)-Orbit.StateEvolution(:,1)); 
            [~, index] = sort(diff);
            index = index(1);

            meas = []; 
            timestamp = [];
            for i = 1:length(Tspan)
                if (dot(Orbit.StateEvolution(index+i,1:3), StateEvolution(i,:)) > 0)
                    meas = [meas; Orbit.StateEvolution(index+i,1:3)];
                    timestamp = [timestamp; Tspan(i)];
                end
            end

            % Apply the probability of detection
            index = logical(randsrc(length(timestamp), 1, [0, 1; 1-obj.PD, obj.PD]));
            timestamp = timestamp(index); 
            meas = meas(index,:); 
            StateEvolution = StateEvolution(index,:);
        end

        % Dynamics 
        function [StateEvolution] = Dynamics(Params, State, Tspan)
            StateEvolution = repmat(State, length(Tspan), 1);
        end

        % Gaussian likelihood function
        function [q] = LikelihoodFunction(Sigma, z, y)
            res = y-z;
            q = exp((0.5*res.'*P^(-1)*res))/sqrt(det(Sigma)*(2*pi)^(size(Sigma,1)));
        end
    end
end