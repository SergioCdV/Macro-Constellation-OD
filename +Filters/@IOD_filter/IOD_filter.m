classdef IOD_filter < Filters.BayesFilter 

    properties
    end

    methods
        % Constructor 
        function [obj] = IOD_filter()
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
        [samples] = UniformTangentQuat(obj, L, m);
    end


end