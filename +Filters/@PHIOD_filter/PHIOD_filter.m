classdef PHIOD_filter < Filters.BayesFilter 

    properties
        % Gravitational problem 
        mu = 3.986e14;                  % Earth's gravitational parameter
        Re = 6378e3;                    % Reference radius of the Earth
        epsilon = -0 * 1.08263e-3;          % J2
        Tc;                             % Characteristic time 

        % Grid properties
        M;
        nu;

        % Resampling
        Jmax = 1e2;  
        PruneThresh = 1e-5;
        MergeThresh = deg2rad(3);
        ResamplingMethod = 'Stratified';

        % State estimation
        N = 1; 
        X;

        % Markov probabilities
        PD = 1; 
        PS = 1;

        % Kalman Filter 
        KF_type = 'UKF-A';
        PD_tol = 1e-6;

        % Phase space delimiters 
        Lmin = 1.03;
        Lmax = 7;
        emax = 0.2;
    end

    methods
        % Constructor 
        function [obj] = PHIOD_filter(myN, myM, myPS, myPD)
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
        end
        
        % Initialization 
        [particles, weights] = Initialization(obj)

        % Bayesian recursion
        [f, X, N, Prior, E] = BayesRecursion(obj, tspan, measurements);
        [PropPrior, sigma] = PropagationStep(obj, last_epoch, new_epoch, Estimator, Prior);
        [Posterior, mean_PD] = CorrectionStep(obj, indices, Measurements, Estimator, PropPrior);
        [planes] = PlanePropagation(obj, planes, step);
    end

    methods 
        % Sampling methods for Gaussians
        [samples] = AffineSampling(obj, m, mu, Sigma);
        [samples] = GibbsSampling(obj, m, mu, Sigma, a); 

        % Clustering and state estimation
        [State] = ParticleState(obj, SensorModality, particle, nu);
        [X] = StateEstimation(obj, N, particles, weights);
        [c, Sigma, index] = QuatClustering(obj, weights, samples);
    end
end