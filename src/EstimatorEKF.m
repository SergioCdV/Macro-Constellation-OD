%% Constellation macro-orbit determination 
% Date: 19/01/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of several EKF Estimators 
% This script provides the function implementing the class EKF

classdef EstimatorEKF
    % Properties
    properties
        Algorithm = 'EKF'   % EKF version in use
        InitFlag = true;

        % State variables
        Clock = 0
        State 
        Sigma
        Measurements

        % State and covariance
        StateDim
        MeasDim

        Q
        R

        StateModel 
        ObservationModel
    end

    % Initialization methods
    methods 
        % State process
        function [obj] = AssignStateProcess(obj, myStateDim, myStateModel)
            obj.StateDim = myStateDim; 
            obj.StateModel = myStateModel;

            if (isempty(obj.R))
                obj.Q = zeros(myStateDim);
            end
        end

        % Observation process assigment 
        function [obj] = AssignObservationProcess(obj, myMeasDim, myObservationProcess)
            obj.MeasDim = myMeasDim;
            obj.ObservationModel = myObservationProcess;

            if (isempty(obj.R))
                obj.R = zeros(myMeasDim);
            end
        end

        % Additive covariances
        function [obj] = AditiveCovariances(obj, Q, R)
            obj.Q = Q;
            obj.R = R;
        end

        % Initial conditions
        function [obj] = InitConditions(obj, State, Sigma)
            obj.State = State;
            obj.Sigma = Sigma;
        end
    end

    % Internal methods
    methods 
        % UKF estimation 
        function [obj] = estimate(obj, time_step, z)
            % Prediction step 
            [state, sigma, y, H, K] = EKF_prediction(obj, time_step);

            % Correction step 
            [state, sigma] = EKF_correction(obj, state, sigma, y, z, H, K); 

            % Clock update 
            obj.Clock = obj.Clock + time_step;

            % Estimator update 
            obj.State = state; 
            obj.Sigma = sigma;
            obj.Measurements = y;
        end

        % UKF prediction 
        function [State, Sigma, y, H, P, K] = EKF_prediction(obj, time_step)
            % Propagation of sigma points 
            [State, Sigma] = obj.StateModel(obj.Q, time_step, obj.State, obj.Sigma);
                
            % Measurement prediction 
            [y, H] = obj.ObservationModel(State);

            % Kalman gain 
            P = (obj.R+H*Sigma*H.');
            K = Sigma*H.'*P^(-1);
        end
        
        % EKF correction
        function [State, Sigma] = EKF_correction(obj, State, Sigma, y, z, H, K)                    
            % Update
            State = State+K*(z-y);
            Sigma = (eye(obj.StateDim)-K*H)*Sigma*(eye(obj.StateDim)-K*H).'+K*obj.R*K.';
        end
    end
end

%% Auxiliary functions
