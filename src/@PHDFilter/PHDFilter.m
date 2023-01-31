%% Constellation macro-orbit determination 
% Date: 19/01/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of the Gaussian Filter Density Hypothesis Filter

classdef PHDFilter
    % Fundamental
    properties
        J = 20                 % Number of initial Gaussian mixtures components
        Jmax = 100             % Maximum number of components
        Jbirth = 0             % Number of birth source components
    
        Mean                   % Mean of the Gaussian components
        Sigma                  % Standard deviation of the Gaussian components
        Weights                % Weights of the mixture

        BirthMeans = 0         % Means of the birth components 
        BirthSigma = 0         % Standard deviation of the Gaussian birth components
        BirthWeights = 0       % Weights of the sources

        PS                     % Probability of surviving
        PD                     % Probability of detection

        PruneThresh = 1e-5     % Threshold to prune components
        MergeThresh = 4        % Threshold to merge two components

        ClutterDensity = 0     % Density of false measurements along the space
        ClutterRate = 0        % Rate of generation of false measurements

        LikelihoodFunction = 1 % Likelihood function to be used
        Domain = 0

        Gaussian;              % Gaussian function to be used
    end

    % Public methods
    methods 
        % General constructor
        function [obj] = PHDFilter(J, Jmax, Mean, Sigma, PS, PD)
            % Principal components
            obj.J = J; 
            obj.Jmax = Jmax; 
            obj.Mean = Mean; 
            obj.Sigma = Sigma; 
            
            % Maximum number components
            obj.Jmax = Jmax; 
        
            % Survivial and detection models 
            obj.PS = PS; 
            obj.PD = PD; 
        end

        % Domain constructor
        function [obj] = DefineDomain(obj, a, b, N)
            obj.Domain = linspace(a,b,N);
        end

        % Define Gaussian functions 
        function [obj] = AddGaussian(obj,myGaussian)
            switch (myGaussian)
                case 'Normal'
                    obj.Gaussian = @(L, sigma, mu)obj.normal(L, sigma, mu);
                case 'Wrapped'
                    obj.Gaussian = @(L, sigma, mu)obj.wrapped_normal(1E-7, L, sigma, mu);
                otherwise
                    error('No implemented Gaussian distribution was selected.')
            end
        end

        % Birth configuration 
        function [obj] = birth(obj, myJbirth, myBirthMeans, myBirthSigma)
            obj.Jbirth = myJbirth; 
            obj.BirthMeans = myBirthMeans;
            obj.BirthSigma = myBirthSigma;

            if (size(myBirthMeans) ~= size(myBirthSigma))
                error('Dimensions mismatch');
            end
        end

        % Pruning heuristic configuration 
        function [obj] = DefinePruning(obj, myPruneThresh, MergeThresh)
            obj.PruneThresh = myPruneThresh; 
            obj.MergeThresh = MergeThresh; 
        end

        % Clutter configuration 
        function [obj] = clutterer(obj, myClutterDensity, myClutterRate)
            obj.ClutterDensity = myClutterDensity;
            obj.ClutterRate = myClutterRate;
        end

        % Likelihood function 
        function [obj] = AssignLikelihood(obj, myLikelihoodFunction)
            obj.LikelihoodFunction = myLikelihoodFunction;
        end
    end

    % Private methods
    methods (Access = protected, Hidden = true)
         % Safety checks
         function [obj] = safety_checks(obj) 
            if (size(obj.Mean,1) == 1 && size(obj.Mean,2) ~= 1)
                warning('Dimensions for the mixture are not consistent.')
                obj.Mean = obj.Mean.';
            end
        
            if (size(obj.Sigma,1) == 1 && size(obj.Sigma,2) ~= 1)
                warning('Dimensions for the mixture are not consistent.')
                obj.Mean = obj.Sigma.';
            end
        
            if (size(obj.BirthMeans,1) == 1 && size(obj.BirthMeans,2) ~= 1)
                warning('Dimensions for the mixture are not consistent.')
                obj.BirthMeans = obj.BirthMeans.';
            end
        
            if (size(obj.BirthSigma,1) == 1 && size(obj.BirthSigma,2) ~= 1)
                warning('Dimensions for the mixture are not consistent.')
                obj.BirthSigma = obj.BirthSigma.';
            end
         end

         % Compute the wrapped normal distribution
         function [f] = normal(obj, L, sigma, mu)
            f = exp(-0.5*(L-mu).^2/sigma);
         end

         % Compute the wrapped normal distribution
         function [f] = wrapped_normal(obj, error_tol, L, sigma, mu)
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

         % Pruning 
         function [w, m, sigma, J] = pruner(obj, j, J, w, m, sigma, i)
            Set = w(i,:) > obj.PruneThresh;
            if (any(Set))
                l = 0;
                while (any(Set))
                    l = l+1; 
                    [~, index] = sort(w(i,Set));
                    max = m(Set,i);
                    P = sigma(Set,i);
                    Mergeable = (max-max(index(end))).^2./sigma(Set,i) <= obj.MergeThresh;
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
            if (size(w(i,:),2) > obj.Jmax)
                [~,index] = sort(w(i,:)); 
                index = index(end-obj.Jmax+1:end);
                w = w(:,index);
                m = m(index,:);
                sigma = sigma(index,:);
                J = obj.Jmax; 
            else
                w = w(:,1:l);
                m = m(1:l,:);
                sigma = sigma(1:l,:);
                J = l; 
            end
         end

         % State estimation
         function [X] = state_estimation(obj, J, N, w, m)
            aux = [];
            for l = 1:J
                if (w(l) > 0.5)
                    aux = [aux; m(l)];
                end
            end
        
            if (N > 0)
                [C, S] = obj.kp_means(N,aux);
                X = [C.'; S.'];
            else
                X = [];
            end
         end

         % K-means 
         function [M, S] = kp_means(obj, N, X)
            % Compute the centroids 
            if (N > size(X,1))
                N = size(X,1);
            end

            [index, M] = kmeans(X,N);
        
            % Compute the covariances
            S = zeros(N,1);
            for i = 1:N
                sigma = (X(index(index == i))-M(i)).^2;
                S(i) = sum(sigma)/sum(index(index == i));
            end
         end
    end
end