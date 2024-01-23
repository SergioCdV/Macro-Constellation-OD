%% Constellation macro-orbit determination 
% Date: 19/01/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of several EKF Estimators 
% This script provides the function implementing the class EKF

classdef EKF < Filters.BayesFilter
    % Properties
    properties
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
        function [obj] = AdditiveCovariances(obj, Q, R)
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
        % EKF estimation 
        function [obj] = BayesStep(obj, time_step, z)
            % Prediction step 
            [state, sigma, y, H, K] = PropagationStep(obj, time_step);

            % Correction step 
            [state, sigma] = CorrectionStep(obj, state, sigma, z); 

            % Clock update 
            obj.Clock = obj.Clock + time_step;

            % Estimator update 
            obj.State = state; 
            obj.Sigma = sigma;
            obj.Measurements = y;
        end

        % EKF prediction 
        [State, Sigma, y, H, P, K] = PropagationStep(obj, time_step);
        
        % EKF correction
        [State, Sigma] = CorrectionStep(obj, State, Sigma, z) ;
    end
end

