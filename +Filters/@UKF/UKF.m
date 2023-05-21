%% Constellation macro-orbit determination 
% Date: 19/01/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of several UKF Estimators 
% This script provides the function implementing the class UKF
% implementation, inlcuding UKF additive and UKF square 

classdef UKF < Filters.BayesFilter
    % Properties
    properties
        % UKF hyperparameters
        Algorithm = 'UKF-A'   % UKF version in use
        beta    
        alpha
        k
        
        Clock = 0
        InitFlag = true
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

    properties (Access = private)
        % Hyperparameters
        W           % UKF weights

        L           % Number of sigma points 
        lambda      % Scaling parameter
        c           % Scaling factor
        sqc         % Square of the scaling factor
    end

    % Initialization methods
    methods 
        % Constructor 
        function [obj] = UKF(Algorithm, beta, alpha, k)
            % Default initialization
            obj.Algorithm = Algorithm; 
            obj.beta = beta; 
            obj.alpha = alpha; 
            obj.k = k; 
        end

        % State process
        function [obj] = AssignStateProcess(obj, myStateDim, myStateModel)
            obj.StateDim = myStateDim; 
            obj.StateModel = myStateModel;

            if (isempty(obj.Q))
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
            switch (obj.Algorithm)
                case 'UKF-S'
                    obj.Q = sqrt(Q);
                    obj.R = sqrt(R);
                otherwise
                    obj.Q = Q;
                    obj.R = R;
            end
        end

        % Initial conditions
        function [obj] = InitConditions(obj, State, Sigma)
            obj.State = State;

            switch (obj.Algorithm)
                case 'UKF-S'
                    if (obj.InitFlag)
                        obj.Sigma = chol(Sigma);
                    end
                    obj.Sigma = chol(Sigma);
                otherwise
                    obj.Sigma = Sigma;
            end

            if (obj.InitFlag)
                obj.InitFlag = false;
            end
        end

        % Initialization 
        function [obj] = Init(obj)
            % Parameters initialization
            obj.L = obj.StateDim;
            obj.lambda = obj.alpha^2*(obj.L+obj.k)-obj.L;
            obj.c = obj.L+obj.lambda;
            obj.sqc = sqrt(obj.c);
            obj.W = [obj.lambda/obj.c*ones(2,1) 0.5/obj.c*ones(2,2*obj.L)]; 
            obj.W(2,1) = obj.W(2,1)+1-obj.alpha^2+obj.beta;
        end

        % UKF estimation 
        function [obj] = estimate(obj, time_step, measurement)
            % Prediction step 
            [sig, state, sigma] = obj.PropagationStep(time_step);

            % Correction step 
            [state, sigma, ~, y] = obj.CorrectionStep(sig, state, sigma, measurement); 

            % Clock update 
            obj.Clock = obj.Clock + time_step;

            % Estimator update 
            obj.State = state; 
            obj.Sigma = sigma;
            obj.Measurements = y;
        end

        % UKF prediction 
        [sigma, State, Sigma] = PropagationStep(obj, time_step);
        
        % UKF correction 
        [State, Sigma, Pmeas, Y] = CorrectionStep(obj, sigma, State, Sigma, z);

        % Sigma points generation
        function [sigma] = sigma_points(obj, State, Sigma)
            switch (obj.Algorithm)
                case 'UKF-S'
                    A = obj.sqc*Sigma.';
                otherwise
                    S = chol(Sigma).'; 
                    A = obj.sqc*S;
            end
            sigma = [State State+A State-A];
        end
    end

    % Private methods 
    methods (Access = private)
        % Measurements reconstruction
        function [Y] = measurements_prediction(obj, y)
            Y = sum(obj.W(1,:).*y,2);
        end

        % General UKF prediction 
        [X, P] = UKFA_prediction(obj, sigma)

        % UKF-A correction
        [State, Sigma, Pyy] = UKFA_correction(obj, sigma, State, Sigma, y, Y, z);

        % UKF-S prediction
        [X, P] = UKFS_prediction(obj, sigma);
        
        % UKF-S correction
        [State, Sigma, Sy] = UKFS_correction(obj, sigma, State, Sigma, y, Y, z);
        
    end
end
