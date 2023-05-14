classdef IOD_filter < Filters.BayesFilter 

    properties
        % Search limits 
        search_limit = [1 7; 0.9 6.9; -6.9 6.9]; 

        % Grid properties
        M; 
        Grid;

        % Tangent grid properties
        L; 

        % Resampling
        Jmax = 3e3; 
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
                obj.Grid = obj.UniformSphere(myM);
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
        end
        
        % Initialization 
        [particles, weights] = Initialization(obj)

        % Bayesian recursion
        [f, X, N] = BayesRecursion(obj, tspan, measurements);
        [PropPrior] = PropagationStep(obj, last_epoch, new_epoch, Prior);
        [Posterior] = CorrectionStep(obj, Measurements, PropPrior, indices);
    end

    methods 
        % Sampling methods for Gaussians
        [samples] = AffineSampling(obj, m, mu, Sigma);
        [samples] = GibbsSampling(obj, m, mu, Sigma, a); 

        % Sampling methods for the sphere
        [samples] = UniformQuat(obj, m); 
        [samples] = UniformSphere(obj, m);

        % Grid
        [born_particles] = Birth(obj);
        [samples] = UniformTangentQuat(obj, L, m, mode);
        [grid] = TransportGrid(samples, mode, post_mode);

        % Clustering and state estimation
        [q] = SteepestQuat(obj, q0, w, dq);
        [c, Sigma, index] = QuatClustering(obj, samples, N); 
        [c, Sigma] = ActionClustering(obj, samples, N);

        % Resampling
        [particles, weights] = Pruning(obj, particles, weights);
        [particles, weights] = Resampling(obj, particles, weights, N);
    end


end