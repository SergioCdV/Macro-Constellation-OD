%% Constellation macro-orbit determination 
% Date: 17/01/2023
% Author: Sergio Cuevas del Valle

%% Kinematic estimator 
% This script provides the function implementing the kinematic estimator
% for the spacecraft intensity function for each orbital plane 

% Inputs: - vector t, the time vector along which the intensity function
%           has to be propagated
%         - array observations, array mx3, which each row is an observation
%         - structure Mixture, containing the definition of the Gaussian
%           Mixture model
%         - class Estimator, the estimator core to be used within the
%           problem

% Outputs:  - 

% Bayesian estimation based on EKF-PHD
function [f, X, N] = kinematic_estimator(obj, t, observations, Estimator)
    % Constants 
    Jb = obj.Jbirth;                   % Number of birth sources
    J = obj.J;                         % Number of mixture components

    mu_b = obj.BirthMeans;             % Mean location of births
    sigma_b = obj.BirthSigma;          % Variance of births location 
    w_b = obj.Birth.Weights;           % Weights of birth mixture

    % False measurements
    Pc = obj.ClutterRate;              % Probability of generating a false measurement
    Vc = obj.ClutterDensity;           % Number of false measurements per orbit
    kappa = Vc*Pc;                     % Clutter distribution

    % Problem estimation parameters
    pd = obj.PD;                       % Detection probability
    ps = obj.PS;                       % Surviving probability

    L = obj.Domain;                    % Estimation domain
    
    % Preallocation 
    m = zeros(J,length(t));             % Mean of the mixtures  
    sigma = zeros(J,length(t));         % Variance of the mixtures
    w = (1/J)*ones(length(t),J);        % Weights in the mixture

    f = zeros(length(L),length(t));     % Intensity function
    M = zeros(length(L),length(t));     % Gaussian functions
    N = zeros(1,length(t));             % Number of targets
    X = cell(length(t),1);              % Multi-target state

    % Safety checks
    obj = safety_checks(obj);

    % Initial state 
    m(:,1) = obj.Mean;                  % Gaussian means    
    sigma(:,1) = obj.Sigma;             % Variance means
    N(1) = 0;                           % Expected number of target

    dim = Estimator.MeasDim;

    % Bayesian estimation
    for i = 1:length(t)
        % Check for measurements
        index = t(i) == observations(:,1);
        j = sum(index);

        % Measurement update
        if (j ~= 0)
            % Collect all measurements at time t_i
            meas = observations(index,:);

            % Preallocation 
            z = zeros(Estimator.MeasDim,J);                                 % Measurements

            switch (Estimator.Algorithm)
                case 'EKF'
                    H = zeros(Estimator.MeasDim,J);                                 % Observation matrix
                    P = zeros(Estimator.MeasDim,Estimator.MeasDim*J);               % Covariance matrix
                    K = zeros(Estimator.StateDim,Estimator.MeasDim*J);              % Kalman gain

                otherwise
                    S = zeros(J, 1, 2*Estimator.StateDim+1);
                    Y = zeros(J, Estimator.MeasDim, 2*Estimator.StateDim+1);
            end

            % Birth proposal 
            J = J+Jb;                       % Current number of Gaussian components
               
            % Step ahead the prior
            if (Jb > 0)
                if (i ~= 1)
                    % Birth
                    m(1:J,i) = [mu_b; m(1:end,i-1)];
                    sigma(1:J,i) = [sigma_b; sigma(1:end,i-1)];
                    w(i,1:J) = [w_b w(i-1,1:end)];   
    
                    % Weight update 
                    w(i,Jb+1:J) = ps*w(i,Jb+1:J);     % Existing weights
                else
                    % Birth
                    m(1:J,i) = [mu_b; m(1:end,i)];
                    sigma(1:J,i) = [sigma_b; sigma(1:end,i)];
                    w(i,1:J) = [w_b w(i,1:end)];    
                end
            end

            for k = 1:Jb
                Dt = 0;
                Estimator.InitFlag = true;
                Estimator = Estimator.InitConditions(m(k,i), sigma(k,i));

                % Prediction
                switch (Estimator.Algorithm)
                    case 'EKF'
                        [m(k,i), sigma(k,i), z(:,k), H(:,k), P(:,1+dim*(k-1):dim*k), K(:,1+dim*(k-1):dim*k)] = Estimator.EKF_prediction(Dt);
                    otherwise
                        [m(k,i), sigma(k,i), z(:,k), S(k,:,:), Y(k,:,:)] = Estimator.UKF_prediction(Dt);
                end
            end

            for k = Jb+1:J
                if (i ~= 1)
                    Dt = t(i)-t(i-1);
                    Estimator = Estimator.InitConditions(m(k,i-1), sigma(k,i-1));
                else
                    Dt = 0;
                    Estimator.InitFlag = true; 
                    Estimator = Estimator.InitConditions(m(k,i), sigma(k,i));
                end

                % Prediction
                switch (Estimator.Algorithm)
                    case 'EKF'
                        [m(k,i), sigma(k,i), z(:,k), H(:,k), P(:,1+dim*(k-1):dim*k), K(:,1+dim*(k-1):dim*k)] = Estimator.EKF_prediction(Dt);
                    otherwise
                        [m(k,i), sigma(k,i), z(:,k), S(k,:,:), Y(k,:,:)] = Estimator.UKF_prediction(Dt);
                end
            end

            % Detections 
            w(i,:) = (1-pd)*w(i,:);

            for k = 1:j
                for l = 1:J
                    % State update with Joseph form for covariance update
                    switch (Estimator.Algorithm)
                        case 'EKF'
                            [m(k*J+l,i), sigma(k*J+l,i)] = Estimator.EKF_correction(m(l,i), sigma(l,i), z(:,l), meas(k,2:end).', H(:,l), K(:,1+dim*(l-1):dim*l));
                            q = likelihood_function(meas(k,2:end).', z(:,l), P(:,1+dim*(l-1):dim*l));
                            
                        otherwise
                            [m(k*J+l,i), sigma(k*J+l,i), P] = Estimator.UKF_correction(m(l,i), sigma(l,i), z(:,l), shiftdim(S(l,:,:)).', shiftdim(Y(l,:,:)), meas(k,2:end).');
                            q = likelihood_function(meas(k,2:end).', z(:,l), P);
                    end
                    
                    if (i ~= 1)
                        w(i,k*J+l) = pd*ps*w(i-1,l)*q;
                    else
                        w(i,k*J+l) = pd*ps*w(i,l)*q;
                    end
                end

                eta = kappa+sum(w(i,1+k*J:(k+1)*J),2);
                if (eta ~= 0)
                    w(i,1+k*J:(k+1)*J) = w(i,1+k*J:(k+1)*J)./eta;
                end
            end

            % Cardinality estimation 
            if (i ~= 1)
                N(i) = round(sum(w(i,J+1:end),2)+(1-pd)*ps*N(i-1));
            else
                N(i) = round(sum(w(i,J+1:end),2)+(1-pd)*ps*N(1));
            end

            % Pruning 
            [w, m, sigma, J] = pruner(obj, j, J, w, m, sigma, i);
        else
            if (i ~= 1)
                % Birth proposal 
                m(1:Jb,i) = m(1:Jb,i-1);
                sigma(1:Jb,i) = sigma(1:Jb,i-1);

                % Number of object
                N(i) = N(i-1);

                Dt = t(i)-t(i-1);
                for k = Jb+1:J
                    Estimator = Estimator.InitConditions(m(k,i-1), sigma(k,i-1));
                    switch (Estimator.Algorithm)
                        case 'EKF'
                            [m(k,i), sigma(k,i), ~, ~, ~] = Estimator.EKF_prediction(Dt);
                        otherwise
                            [m(k,i), sigma(k,i), ~, ~, ~] = Estimator.UKF_prediction(Dt);
                    end
                end
            end
        end

        % Intensity function at time t_i
        m(:,i) = mod(m(:,i), 2*pi);
        for k = 1:J
            M(:,k) = obj.wrapped_normal(L.', sigma(k,i), m(k,i), 1e-7);
        end
        f(:,i) = sum(w(i,:).*M(:,1:size(w,2)),2);

        % Multi-target state
        aux = [];
        for l = 1:J
            if (w(i,l) > 0.5)
                aux = [aux; m(l,i)];
            end
        end

        if (N(i) > 0)
            [C, S] = obj.kp_means(N(i),aux);
            X{i} = [C.'; S.'];
        else
            X{i} = [];
        end
    end
end
