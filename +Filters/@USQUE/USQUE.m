%% Constellation macro-orbit determination 
% Date: 19/01/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of several USQUE Estimator 
% This script provides the function implementing the class USQQUE
% implementation, inlcuding USQUE additive and UKF square 

classdef USQUE < Filters.BayesFilter
    % Properties
    properties
        % UKF hyperparameters
        Algorithm = 'UKF-A'   % UKF version in use
        beta = 2;  
        alpha
        k = 0;
        W                     % UKF weights
        
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

        a 
        f
    end

    properties (Access = private)
        % Hyperparameters
        L           % Number of sigma points 
        lambda      % Scaling parameter
        c           % Scaling factor
        sqc         % Square of the scaling factor
    end

    % Initialization methods
    methods 
        % Constructor 
        function [obj] = USQUE(Algorithm, beta, alpha, k, a, f)
            % Default initialization
            obj.alpha = alpha; 

            if (exist('Algorithm', 'var'))
                obj.Algorithm = Algorithm; 
            end

            if (exist('beta', 'var'))
                obj.beta = beta; 
            end

            if (exist('k', 'var'))
                obj.k = k; 
            end 

            if (exist('a', 'var'))
                obj.a = a; 
            else
                obj.a = 1;
            end 

            if (exist('f', 'var'))
                obj.f = f; 
            else
                obj.f = 2 * (obj.a + 1); 
            end 
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
                        Sigma = chol(Sigma).';
                        obj.InitFlag = false;
                    end
            end

            obj.Sigma = Sigma;
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
            obj.State = [zeros(3,1); state(4:6,:)]; 
            obj.Sigma = sigma;
            obj.Measurements = y;
        end

        % UKF prediction 
        [sigma, State, Sigma] = PropagationStep(obj, time_step);
        
        % UKF correction 
        [State, Sigma, Pmeas, Y] = CorrectionStep(obj, sigma, State, Sigma, z);

        % Sigma points generation
        function [sigma] = sigma_points(obj, State, Sigma, Q)
            Sigma = Sigma + Q;

            switch (obj.Algorithm)
                case 'UKF-S'
                    A = Sigma.';
                otherwise
                    A = chol(Sigma).'; 
            end

            A = obj.sqc*A;
            sigma = [State State+A State-A];
        end
    end

    % Private methods 
    methods (Access = private)
        % Measurements reconstruction
        function [Y] = measurements_prediction(obj, y)
            Y = sum(obj.W(1,:).*y,2);
        end

        % UKF-A prediction 
        [X, P] = UKFA_prediction(obj, sigma)

        % UKF-A correction
        [State, Sigma, Pyy] = UKFA_correction(obj, sigma, State, Sigma, y, Y, z);        
    end
end
