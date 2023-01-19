%% Constellation macro-orbit determination 
% Date: 19/01/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of several UKF Estimators 
% This script provides the function implementing the class UKF
% implementation, inlcuding UKF additive and UKF square 

classdef EstimatorUKF
    % Properties
    properties
        % UKF hyperparameters
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
        Algorithm = 'UKF-A'   % UKF version in use
        W                     % UKF weights

        L           % Number of sigma points 
        lambda      % Scaling parameter
        c           % Scaling factor
        sqc         % Square of the scaling factor
    end

    % Initialization methods
    methods 
        % Constructor 
        function [obj] = EstimatorUKF(Algorithm, beta, alpha, k)
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
        function [obj] = AditiveCovariances(obj, Q, R)
            obj.Q = Q;
            obj.R = R;
        end

        % Initial conditions
        function [obj] = InitConditions(obj, State, Sigma)
            obj.State = State;

            switch (obj.Algorithm)
                case 'UKF-S'
                    if (obj.InitFlag)
                        obj.Sigma = chol(Sigma).';
                    end
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
    end

    % Internal methods
    methods 
        % UKF estimation 
        function [obj] = estimate(obj, time_step, measurement)
            % Prediction step 
            [state, sigma, y, Y, sig] = UKF_prediction(obj, time_step);

            % Correction step 
            [state, sigma, ~] = UKF_correction(obj, state, sigma, y, sig, Y, measurement); 

            % Clock update 
            obj.Clock = obj.Clock + time_step;

            % Estimator update 
            obj.State = state; 
            obj.Sigma = sigma;
            obj.Measurements = y;
        end

        % UKF prediction 
        function [State, Sigma, Measurements, sigma, y] = UKF_prediction(obj, time_step)
            % Generate sigma points 
            sigma = sigma_points(obj);

            % Propagation of sigma points 
            sigma = obj.StateModel(time_step, sigma);
    
            % State and covariance prediction 
            switch (obj.Algorithm)
                case 'UKF-A'
                    [State, Sigma] = UKFA_prediction(obj, sigma);
                case 'UKF-S'
                    [State, Sigma] = UKFS_prediction(obj, sigma);
            end

            sigma = sigma_points(obj);
            
            % Measurement prediction 
            y = obj.ObservationModel(sigma);
            Measurements = measurements_prediction(obj, y);
        end
        
        % UKF correction 
        function [State, Sigma, Pmeas] = UKF_correction(obj, State, Sigma, y, sigma, Y, z)
            % State and covariance prediction 
            switch (obj.Algorithm)
                case 'UKF-A'
                    [State, Sigma, Pmeas] = UKFA_correction(obj, State, Sigma, y, sigma, Y, z);
                case 'UKF-S'
                    [State, Sigma, Sy] = UKFS_correction(obj, State, Sigma, y, sigma, Y, z);
                    Pmeas = Sy.'*Sy;
            end
        end
    end

    % Private methods 
    methods (Access = private)
        % Sigma points generation
        function [sigma] = sigma_points(obj)
            switch (obj.Algorithm)
                case 'UKF-S'
                    A = obj.sqc*obj.Sigma;
                otherwise
                    A = obj.sqc*sqrt(obj.Sigma);
            end
            sigma = [obj.State obj.State+A obj.State-A];
        end

        % Measurements reconstruction
        function [Y] = measurements_prediction(obj, y)
            Y = sum(obj.W(1,:).*y,2);
        end

        % General UKF prediction 
        function [X, P] = UKFA_prediction(obj, sigma)
            % State prediction
            X = sum(obj.W(1,:).*sigma,2);
        
            % Covariance prediction
            P = zeros(obj.StateDim); 
            for i = 1:size(obj.W,2)
                P = P+obj.W(2,i).*(sigma(:,i)-X).*(sigma(:,i)-X).';
            end
            P = P+obj.Q;
        end

        % UKF-S prediction
        function [X, P] = UKFS_prediction(obj, sigma)
            % State prediction
            X = sum(obj.W(1,:).*sigma,2);

            % Covariance prediction
            P = [sqrt(obj.W(2,2))*(sigma(2:end)-X) sqrt(obj.Q)];
            [~, P] = qr(P.',0);
            P = cholupdate(P, obj.W(2,1)*(sigma(:,1)-X));
        end

        % UKF-A correction
        function [State, Sigma, Pyy] = UKFA_correction(obj, State, Sigma, y, sigma, Y, z)
            Pyy = 0; 
            for i = 1:size(obj.W,2)
                Pyy = Pyy+obj.W(2,i)*(Y(:,i)-y)*(Y(:,i)-y).';
            end
            Pyy = Pyy+obj.R;
        
            Pxy = 0; 
            for i = 1:size(obj.W,2)
                Pxy = Pxy+obj.W(2,i)*(sigma(:,i)-State)*(Y(:,i)-y).';
            end
            
            % Kalman gain
            K = Pxy*(Pyy^(-1));
        
            % Update
            State = State+K*(z-y);
            Sigma = Sigma-K*Pyy*K.';
        end
        
        % UKF-S correction
        function [State, Sigma, Sy] = UKFS_correction(obj, State, Sigma, y, sigma, Y, z)
            % Covariance computation
            Sy = [sqrt(obj.W(2,2))*(Y(2:end)-y) sqrt(obj.R)];
            [~, Sy] = qr(Sy.', 0);
            Sy = cholupdate(Sy, obj.W(2,1)*(Y(:,1)-y));
                 
            Pxy = 0;
            for i = 1:size(obj.W,2)
                Pxy = Pxy+obj.W(2,i).*(sigma(:,i)-State)*(Y(:,i)-y).';
            end

            % Kalman update
            K = Pxy/Sy/Sy.';
        
            % Update
            State = State+K*(z-y);
            U = K*Sy.';
            Sigma = cholupdate(Sigma, U(:,1), '-');            
        end
    end
end

%% Auxiliary functions
