classdef (Abstract) BayesFilter 

    properties
    end

    methods
        % Constructor 
        function [obj] = BayesFilter()
        end
        
        % Bayesian recursion
        [prior_m] = propagation_step(prior);
        [posterior] = correction_step(measurements, prior_m);
    end
end