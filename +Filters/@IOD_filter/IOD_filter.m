classdef IOD_filter < Filters.BayesFilter 

    properties
        % Gravitational problem 
        mu = 3.986e14;                  % Earth's gravitational parameter
        Re = 6378e3;                    % Reference radius of the Earth
        epsilon = -1.08263e-3;          % J2
        Tc;                             % Characteristic time 

        % Search limits 
        search_limit = [1 7; 0.9 6.9; -6.9 6.9]; 

        % Grid properties
        M; 
        Grid;

        % Tangent grid properties
        L; 

        % Resampling
        Jmax = 1e4; 
        RevThresh = 1e-4; 
        PruneThresh = 1e-5;
        MergeThresh = deg2rad(3);
        ResamplingMethod = 'Systematic';

        % State estimation
        N = 10; 
        X;

        % Markov probabilities
        PD = 1; 
        PS = 1;

        % Covariance regularization
        PD_tol = 1e-6;
    end

    properties (Access = private)
        W; 
        nu; 
        dt;
    end

    methods
        % Constructor 
        function [obj] = IOD_filter(myM, myL, myN, myPD, myPS)
            % Contruct the filter sample space dimensions
            if (exist('myL', 'var'))
                obj.L = myL;
            end

            if (exist('myM', 'var'))
                obj.M = myM;
    
                % Construct the grid 
                obj.Grid = QuaternionAlgebra.UniformSphere(myM);
            end

            if (exist('myN', 'var'))
                obj.N = myN;
            end

%             if (exist('myNs', 'var'))
%                 obj.Ns = 10;
%             end

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
        [f, X, N, E] = BayesRecursion(obj, tspan, measurements);
        [PropPrior] = PropagationStep(obj, last_epoch, new_epoch, Prior);
        [Posterior] = CorrectionStep(obj, indices, Measurements, PropPrior);
    end

    methods 
        % Sampling methods for Gaussians
        [samples] = AffineSampling(obj, m, mu, Sigma);
        [samples] = GibbsSampling(obj, m, mu, Sigma, a); 

        % Grid
        [born_particles] = Birth(obj);
        [samples] = UniformTangentQuat(obj, L, m, mode);
        [grid] = TransportGrid(samples, mode, post_mode);

        % Clustering and state estimation
        [State] = ParticleState(obj, SensorModality, particle, nu);
        [c, Sigma, index] = QuatClustering(obj, weights, samples); 
        [c, Sigma] = ActionClustering(obj, samples, N);

        % Resampling
        [particles, weights] = Pruning(obj, particles, weights);
        [particles, weights] = Resampling(obj, particles, weights, N);
        [samples, a] = PerifocalQuatSampling(obj, particles);
    end

    methods (Access = private)
        [Y] = quadrature(obj, x);
    end


end