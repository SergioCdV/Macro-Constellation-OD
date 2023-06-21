classdef (Abstract) AbstractSensor
    
    properties
        State = [];             % State of the observer
        Measurements = {};      % Measurements acquired by the observer

        InstrumentParams;       % Instrument parameters

        InitialEpoch;           % Initial epoch
        PropagatedEpoch;        % Propgated epoch
        CurrentEpoch;           % Current epoch

        Sigma;                  % Noise covariance
        PD;                     % Probability of detection
        SR = 1;                 % Sampling rate [s]
        PC = 0;                 % Clutter probability
        NC = 0;                 % Average number of false measurements

        StateDim;               % Dimension of the state vector
        MeasDim;                % Dimension of the measurement vector
    end
    
    methods
        function obj = AbstractSensor(myStateDim, myMeasDim, myInitialEpoch, myInitialState, mySigma, myPD, mySR)
            % Dimensionality checks
            obj.StateDim = myStateDim; 
            obj.MeasDim = myMeasDim; 
            
            % Initial state conditions
            if (myInitialEpoch >= 0)
                obj.InitialEpoch = myInitialEpoch; 
            end

            if (size(myInitialState,2) == obj.StateDim)
                obj.State = myInitialState;
            else
                error('No valid observer state has been addded.');
            end

            % Covariances 
            try chol(mySigma);
                if (size(mySigma,1) == 1)
                    obj.Sigma = mySigma * eye(obj.StateDim);
                else
                    obj.Sigma = mySigma;
                end
            catch
                error('Input covariance matrix is not positive definite.'); 
            end

            if (any(size(obj.Sigma) ~= obj.MeasDim))
                error('Input covariance matrix is not appropriately shaped.');
            end

            if (myPD >= 0 && myPD <= 1)
                obj.PD = myPD;
            else
                error('No valid detection probability has been input.'); 
            end

            if (exist('mySR', 'var'))
                if (mySR)
                    obj.SR = mySR;
                end
            end
        end

        [q] = LikelihoodFunction(obj, Sigma, z, y);
        [obj] = ConfigClutter(myPC, myNC);
    end

    methods (Access = private)
        [Tspan, StateEvolution] = Dynamics(obj, Epoch, State, Tspan);
    end
end

