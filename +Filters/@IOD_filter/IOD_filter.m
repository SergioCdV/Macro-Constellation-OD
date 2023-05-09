classdef IOD_filter < Filters.BayesFilter 

    properties
        % Grid properties
        M; 
        Grid;

        % Tangent grid properties
        L; 

        % Modes
        N; 
        modes;

        % Markov probabilities
        PD = 1; 
        PS = 1;
    end

    methods
        % Constructor 
        function [obj] = IOD_filter(myM, myL, myPD, myPS)
            % Contruct the filter sample space dimensions
            if (exist('myL', 'var'))
                obj.L = myL;
            end

            if (exist('myM', 'var'))
                obj.M = myM;
    
                % Construct the grid 
                obj.Grid = obj.UniformSphere(myM);
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
        end
        
        % Bayesian recursion
        [prior_m] = propagation_step(prior);
        [posterior] = correction_step(measurements, prior_m);
    end

    methods 
        % Sampling methods for Gaussians
        [samples] = AffineSampling(obj, m, mu, Sigma);
        [samples] = GibbsSampling(obj, m, mu, Sigma);

        % Sampling methods for the sphere
        [samples] = UniformQuat(obj, m); 
        [samples] = UniformSphere(obj, m);

        % Grid
        [samples] = UniformTangentQuat(obj, L, m, mode);
        [grid] = TransportGrid(samples, mode, post_mode);

        % Clustering and state estimation
        [c, Sigma, index] = QuatClustering(obj, samples, N); 
        [c, Sigma] = ActionClustering(obj, samples, N);
    end


end