%% Constellation macro-orbit determination 
% Date: 06/02/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of an GibbsObserver object

classdef GibbsObserver
    properties
        Type = 'Gibbs';         % Class of observer

        State = [];             % State of the observer
        Measurements = {};      % Measurements acquired by the observer

        InstrumentParams;       % Instrument parameters
        Sigma;                  % Covariance of the instrument

        InitialEpoch;           % Initial epoch
        PropagatedEpoch;        % Propgated epoch
        CurrentEpoch;           % Current epoch

        PD;                     % Probability of detection
    end

    % Private properties
    properties (Access = private)
        StateDim = 3;           % Dimension of the state vector
        MeasDim = 3;            % Dimension of the measurement vector
    end

    methods 
        % Add probability of detection 
        function [obj] = probability_detection(obj, myPD)
            obj.PD = myPD;
        end

        % Add instrument covariance 
        function [obj] = AddCovariance(obj, mySigma)
            try chol(mySigma);
                if (size(mySigma,1) == 1)
                    obj.Sigma = mySigma * eye(obj.StateDim);
                else
                    obj.Sigma = mySigma;
                end
            catch
                error('Input covariance matrix is not positive definite.'); 
            end
        end

        % Add an initial state 
        function [obj] = AddInitialState(obj, myInitialEpoch, myInitialState)
            if (myInitialEpoch > 0)
                obj.InitialEpoch = myInitialEpoch; 
            end

            if (size(myInitialState,2) == obj.StateDim)
                obj.State = myInitialState;
            else
                error('No valid observer state has been addded.');
            end
        end

        % Observations
        function [timestamp, meas, StateEvolution] = Observe(obj, Orbit, Epoch)
            % Propagate the orbit 
            AuxOrbit = Orbit.SetCurrentEpoch(Epoch).Propagate().ChangeStateFormat('Cartesian');

            % Propagate the observer and take the measurements
            AuxOrbit = AuxOrbit.Normalize(false, 1);
            Tspan = AuxOrbit.StateEvolution(:,1);
            [Tspan, StateEvolution] = obj.Dynamics(obj.State, Tspan+AuxOrbit.InitialEpoch); 

            index = Tspan >= obj.InitialEpoch;
            AuxOrbitEvolution = AuxOrbit.StateEvolution(index,2:end);
            Tspan = Tspan(index);

            if (~isempty(Tspan) && ~isempty(StateEvolution))
                % Restrict the orbit state evolution    
                meas = []; 
                timestamp = [];
                
                [timestamp, meas] = obj.ObservationProcess(Tspan, AuxOrbitEvolution, StateEvolution);
    
                % Apply the probability of detection
                index = logical(randsrc(length(timestamp), 1, [0, 1; 1-obj.PD, obj.PD]));
                timestamp = timestamp(index); 
                meas = meas(index,:); 
                StateEvolution = StateEvolution(index,:);
            else
                timestamp = []; 
                meas = []; 
                StateEvolution = [];
            end
        end

        % Observation process 
        function [t, meas] = ObservationProcess(obj, Tspan, Orbit, StateEvolution)
            % Preallocaton
            meas = []; 
            t = []; 

            % Observation
            for i = 1:length(Tspan)
                if (dot(Orbit(i,1:3), StateEvolution(i,:)) > 0)
                    meas = [meas; Orbit(i,1:3)];
                    t = [t; Tspan(i)];
                end
            end
        end

        % Gaussian likelihood function
        function [q] = LikelihoodFunction(obj, Sigma, z, y)
            res = y-z;
            q = exp((0.5*res.'*P^(-1)*res))/sqrt(det(Sigma)*(2*pi)^(size(Sigma,1)));
        end
    end
    
    methods (Access = private)
        % Dynamics 
        function [Tspan, StateEvolution] = Dynamics(obj, State, Tspan)
            % Check if an initial state has been given 
            if (~isempty(obj.State))
                % Check the tspan 
                index = Tspan >= obj.InitialEpoch;
    
                % Copy the state
                StateEvolution = repmat(State, length(Tspan(index)), 1);
            else
                error('No obsever state has been given.');
            end
        end
    end
end