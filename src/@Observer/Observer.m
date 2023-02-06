%% Constellation macro-orbit determination 
% Date: 06/02/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of an Observer object

classdef Observer
    % Properties 
    properties 
        Type;                   % Class of observer

        State;                  % State of the observer
        Measurements = {};      % Measurements acquired by the observer

        InstrumentParams;       % Instrument parameters

        InitialEpoch;           % Initial epoch
        PropagatedEpoch;        % Propgated epoch
        CurrentEpoch;           % Current epoch

        Likelihood;             % Instrument likelihood function

        Dynamics;               % Dynamics of the observer
        ObservationModel;       % Observation function

        Tspan;                  % Considered time span
    end

    % Private properties
    properties (Access = private)
        StateDim;               % Dimension of the state vector
        MeasDim;                % Dimension of the measurement vector
    end

    % Define the observer 
    methods 
        % Constructor
        function [obj] = Observer(myType)
            switch (myType)
                case 'Ideal telescope'
                    obj = IdealTelescope();

                case 'Radec telescope'
                case 'Radar'
                otherwise
                    error('No valid observer type was selected');
            end
        end

        % Set initial epoch
        function [obj] = SetInitialEpoch(obj, myInitialEpoch)
            if (myInitialEpoch > 0)
                obj.InitialEpoch = myInitialEpoch;
            end
        end

        % Set current epoch
        function [obj] = SetCurrentEpoch(obj, myCurrentEpoch)
            if (myCurrentEpoch > 0)
                obj.CurrentEpoch = myCurrentEpoch;
            end
        end

        % Set covariance 
        function [obj] = SetInstrumentParameters(obj, myInstrumentParameters)
            obj.InstrumentParams = myInstrumentParameters;
            obj.Likelihood = @(z, y)obj.Likelihood(obj.Params, z, y);
        end

        % Define the initial state of the observer 
        function [obj] = SetState(obj, myState)
            if (size(myState,1) ~= obj.StateDim)
                error('Dimension mismatch for the proposed state vector.');
            else
                obj.State = myState;
            end
        end
                
        % Observe an orbit or constellation
        function [obj] = Observe(obj, OrbitSet)
            % Observe the orbit set
            if (~isempty(OrbitSet))
                for i = 1:size(OrbiSet,1)
                    [timestamps, meas, StateEvolution] = obj.ObservationModel(OrbitSet{i,:});
                    obj.Measurements = [obj.Measurements; {timestamps, meas, StateEvolution, @(Orbit)obj.ObservationModel(Orbit), @(y)obj.Likelihood(obj.InstrumentParams, meas, y)}];
                end
            end

            % Set the new propagated epoch 
            obj.PropagatedEpoch = obj.CurrentEpoch;
            obj.State = StateEvolution(end,:);
        end
    end

    % Private methods 
    methods (Access = private)
        % Propagate the observer state
        function [StateEvolution] = Propagate(obj, tspan)
            StateEvolution = obj.Dynamics(obj.State, tspan);
        end
    end
end