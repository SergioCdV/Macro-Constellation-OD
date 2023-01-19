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
function [L, f, N, X] = kinematic_estimator_EKF(t, observations, Mixture, Estimator)
    % Constants 
    prune_thresh = Mixture.Th.Prune;    % Threshold to prune components
    merge_thresh = Mixture.Th.Merge;    % Threshold to merge components
    Jmax = Mixture.Jmax;                % Maximum number of components 
    Jb = Mixture.Birth.J;               % Number of birth sources
    J = Mixture.J;                      % Number of mixture components

    mu_b = Mixture.Birth.Mean;          % Mean location of births
    sigma_b = Mixture.Birth.Sigma;      % Variance of births location 
    w_b = Mixture.Birth.Weights;        % Weights of birth mixture

    % False measurements
    Pc = Mixture.Clutter.Rate;          % Probability of generating a false measurement
    Vc = Mixture.Clutter.Density;       % Number of false measurements per orbit
    kappa = Vc*Pc;                      % Clutter distribution

    % Problem estimation parameters
    pd = Mixture.Probabilities(1);          % Detection probability
    ps = Mixture.Probabilities(2);          % Surviving probability
    
    % Angular domain
    Npart = 1e3;                        % Discretization of the angular domain
    L = linspace(0,2*pi,Npart);         % Angular domain draws

    % Preallocation 
    m = zeros(J,length(t));             % Mean of the mixtures  
    sigma = zeros(J,length(t));         % Variance of the mixtures
    w = (1/J)*ones(length(t),J);        % Weights in the mixture

    f = zeros(length(L),length(t));     % Intensity function
    M = zeros(length(L),length(t));     % Gaussian functions
    N = zeros(1,length(t));             % Number of targets
    X = cell(length(t),1);              % Multi-target state

    Mixture = safety_checks(Mixture);

    % Initial state 
    m(:,1) = Mixture.Mean;              % Gaussian means    
    sigma(:,1) = Mixture.Sigma;         % Variance means
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
            H = zeros(Estimator.MeasDim,J);                                 % Observation matrix
            P = zeros(Estimator.MeasDim,Estimator.MeasDim*J);               % Covariance matrix
            K = zeros(Estimator.StateDim,Estimator.MeasDim*J);              % Kalman gain

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
                Estimator = Estimator.InitConditions(m(k,i), sigma(k,i));

                % Prediction
                [m(k,i), sigma(k,i), z(:,k), H(:,k), P(:,1+dim*(k-1):dim*k), K(:,1+dim*(k-1):dim*k)] = Estimator.EKF_prediction(Dt);
            end

            for k = Jb+1:J
                if (i ~= 1)
                    Dt = t(i)-t(i-1);
                    Estimator = Estimator.InitConditions(m(k,i-1), sigma(k,i-1));
                else
                    Dt = 0;
                    Estimator = Estimator.InitConditions(m(k,i), sigma(k,i));
                end

                % Prediction
                [m(k,i), sigma(k,i), z(:,k), H(:,k), P(:,1+dim*(k-1):dim*k), K(:,1+dim*(k-1):dim*k)] = Estimator.EKF_prediction(Dt);
            end

            % Detections 
            w(i,:) = (1-pd)*w(i,:);

            for k = 1:j
                for l = 1:J
                    % State update with Joseph form for covariance update
                    [m(k*J+l,i), sigma(k*J+l,i)] = Estimator.EKF_correction(m(l,i), sigma(l,i), z(:,l), meas(k,2:end).', H(:,l), K(:,1+dim*(l-1):dim*l));

                    % Weight update 
                    q = likelihood_function(meas(k,2:end).', z(:,l), P(:,1+dim*(l-1):dim*l));
                    
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
            Set = w(i,:) > prune_thresh;
            if (any(Set))
                l = 0;
                while (any(Set))
                    l = l+1; 
                    [~, index] = sort(w(i,Set));
                    max = m(Set,i);
                    P = sigma(Set,i);
                    Mergeable = (max-max(index(end))).^2./sigma(Set,i) <= merge_thresh;
                    aux = w(i,Set);
                    w(i,l) = sum(aux(Mergeable));
    
                    m(l,i) = dot(aux(Mergeable),max(Mergeable))/w(i,l);
                    sigma(l,i) = dot(aux(Mergeable),(P(Mergeable)+(m(l,i)-max(Mergeable)).^2))/w(i,l);
                    
                    q = 1;
                    for k = 1:length(Set)
                        if (Set(k))
                            Set(k) = ~Mergeable(q);
                            q = q+1;
                        end
                    end
                end
            else
               l = (j+1)*J;
            end

            % Check if there are more than Jmax components 
            if (size(w(i,:),2) > Jmax)
                [~,index] = sort(w(i,:)); 
                index = index(end-Jmax+1:end);
                w = w(:,index);
                m = m(index,:);
                sigma = sigma(index,:);
                J = Jmax; 
            else
                w = w(:,1:l);
                m = m(1:l,:);
                sigma = sigma(1:l,:);
                J = l; 
            end
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
                    [m(k,i), sigma(k,i), ~, ~, ~] = Estimator.EKF_prediction(Dt);
                end
            end
        end

        % Intensity function at time t_i
        m(:,i) = mod(m(:,i), 2*pi);
        for k = 1:J
            M(:,k) = wrapped_normal(L.', sigma(k,i), m(k,i), 1e-7);
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
            [C, S] = kp_means(N(i),aux);
            X{i} = [C.'; S.'];
        else
            X{i} = [];
        end
    end
end

%% Auxiliary functions 
% Safety checks
function [Mixture] = safety_checks(Mixture)
    % Safety checks 
    if (size(Mixture.Mean,1) == 1 && size(Mixture.Mean,2) ~= 1)
        warning('Dimensions for the mixture are not consistent.')
        Mixture.Mean = Mixture.Mean.';
    end

    if (size(Mixture.Sigma,1) == 1 && size(Mixture.Sigma,2) ~= 1)
        warning('Dimensions for the mixture are not consistent.')
        Mixture.Mean = Mixture.Sigma.';
    end

    if (size(Mixture.Birth.Mean,1) == 1 && size(Mixture.Birth.Mean,2) ~= 1)
        warning('Dimensions for the mixture are not consistent.')
        Mixture.Mean = Mixture.Birth.Mean.';
    end

    if (size(Mixture.Birth.Sigma,1) == 1 && size(Mixture.Birth.Sigma,2) ~= 1)
        warning('Dimensions for the mixture are not consistent.')
        Mixture.Mean = Mixture.Birth.Sigma.';
    end
end

% Compute the wrapped normal distribution
function [f] = wrapped_normal(L, sigma, mu, error_tol)
    % Compute the error bound 
    n(1) = max(1+sqrt(-log(4*pi^3*error_tol^2)*sigma),1+sqrt(sigma/2)/pi);
    n(2) = max(sqrt(-log(2*pi^2*sigma*error_tol^2)/sigma),sqrt(2)/pi);
    n = ceil(n);

    f = zeros(length(L),1);

    % Compute the wrapped normal distributions
    if (min(n) == n(1))
        % Compute the wrapped normal
        for i = -n:n
            f = f+exp(-0.5*(L-mu+2*pi*i).^2/sigma);
        end
        f = f/sqrt(sigma*2*pi);
    else
        rho = exp(-sigma/2);
        for i = 1:n
            f = f+rho^(i^2)*cos(i*(L-mu));
        end
        f = (1+2*f)/(2*pi);
    end
end

% Compute the likelihood function 
function [q] = likelihood_function(z, m, P)
    q = exp(-0.5*(z-m).'*P^(-1)*(z-m))/sqrt(det(P)*(2*pi)^size(P,1));
end

% K-means 
function [M, S] = kp_means(N,X)
    % Compute the centroids 
    [index, M] = kmeans(X,N);

    % Compute the covariances
    S = zeros(N,1);
    for i = 1:N
        sigma = (X(index(index == i))-M(i)).^2;
        S(i) = sum(sigma)/sum(index(index == i));
    end
end
    % Constants 
    tol = 1e-3; 

    % Initialize at random the clusters 
    P = cell(N,1);            % Partitions
    M = zeros(N,1);           % Means
    S = zeros(N,1);           % Covariance
    
    n = 0;
    i = 1; 
    while (n < N)
        if (randi([0 1]))
            n = n+1; 
            M(n) = X(i);
        elseif (i+1 == length(X))
            i = 0;
        end
        i = i+1;
    end
    j = 2; 

    % Compute the means 
    GoOn = true; 
    iter = 1; 
    maxIter = 100;
    while (GoOn && iter < maxIter)
        % Compute the partitions 
        for i = 1:length(X)
            min_index = 1;
            d = 1e7;
            for j = 1:N
                aux = norm(X(i)-M(j));
                if (aux < d)
                    min_index = j; 
                    d = aux;
                end
            end

            P{min_index} = [P{min_index} X(i)];
        end

        % Compute the means
        for i = 1:N
            M(i) = mean(P{i});
        end
        j = j+1;

        % Convergence analysis
        if (abs(cost - past_cost) < tol)
        else
            iter = iter+1;
        end
    end
    GoOn = ~GoOn;

    % Compute the covariances
end