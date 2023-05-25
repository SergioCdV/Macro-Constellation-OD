classdef MTT_filter < Filters.BayesFilter 

    properties
        % Gravitational problem 
        mu = 3.986e14;                  % Earth's gravitational parameter
        Re = 6378e3;                    % Reference radius of the Earth
        epsilon = -1.08263e-3;          % J2
        Tc;                             % Characteristic time 

        % Grid properties
        M;
        nu;

        % Resampling
        Jmax = 1e3;  
        PruneThresh = 1e-5;
        MergeThresh = deg2rad(3);
        ResamplingMethod = 'Systematic';

        % State estimation
        planes;
        N = 1; 
        X;

        % Markov probabilities
        PD = 1; 
        PS = 1;

        % Kalman Filter 
        KF_type = 'UKF-A';
    end

    methods
        % Constructor 
        function [obj] = MTT_filter(myN, myM, myPD, myPS, myPlanes)
            % Contruct the filter sample space dimensions
            if (exist('myM', 'var'))
                obj.M = myM;
            end

            if (exist('myN', 'var'))
                obj.N = myN;
            end

            if (exist('myPD', 'var'))
                if (myPD >= 0 && myPD <= 1)
                    obj.PD = myPD;
                else
                    error('No valid detection probability has been input.'); 
                end
            end

            if (exist('myPS', 'var'))
                if (myPS >= 0 && myPS <= 1)
                    obj.PS = myPS;
                else
                    error('No valid survival probability has been input.'); 
                end
            end

            obj.Tc = sqrt(obj.Re^3/obj.mu);  

            obj.planes = myPlanes;
        end
        
        % Initialization 
        [particles, weights] = Initialization(obj)

        % Bayesian recursion
        [f, X, N, Prior] = BayesRecursion(obj, tspan, measurements);
        [PropPrior, sigma] = PropagationStep(obj, last_epoch, new_epoch, Estimator, Prior);
        [Posterior] = CorrectionStep(obj, indices, Measurements, Estimator, PropPrior);
        [Prior] = PlanePropagation(obj, last_epoch, prop_epoch, planes);
    end

    methods 
        % Wrapped normal 
        [f] = wrapped_normal(obj, error_tol, M, mu, sigma); 

        % Sampling methods for Gaussians
        [samples] = AffineSampling(obj, m, mu, Sigma);
        [samples] = GibbsSampling(obj, m, mu, Sigma, a); 

        % Clustering and state estimation
        [X] = StateEstimation(obj, N, particles, weights);

        % Resampling
        [particles, weights] = Pruning(obj, particles, weights);
    end

end