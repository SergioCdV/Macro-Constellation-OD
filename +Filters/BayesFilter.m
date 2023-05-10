classdef (Abstract) BayesFilter 

    properties
    end

    methods
        % Constructor 
        function [obj] = BayesFilter()
        end

        % Initialization 
        [obj] = initialization(obj);
        
        % Bayesian recursion
        [PropPrior] = PropagationStep(obj, last_epoch, new_epoch, Prior);
        [Posterior] = CorrectionStep(obj, Measurements, PropPrior);
    end
end