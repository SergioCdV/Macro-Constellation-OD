

function [particles, weights] = Initialization(obj)
    % Sanity checks 
    if (isempty(obj.N) || isempty(obj.L) || isempty(obj.M))
        error('The internal filter configuration has not been completed.');
    else  
        % Anomaly distribution initialization 
        mu = linspace(0, 2*pi, J);
        sigma = 1e-3 * ones(1, J);

        % Particles as wrapped Gaussian kernels 
        particles = [mu; sigma];

        % Uniform generation of the weights 
        weights = repmat(obj.N/size(mu,2), 1, size(mu,2));
    end
end