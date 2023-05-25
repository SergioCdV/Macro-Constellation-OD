

function [particles, weights] = Initialization(obj)
    % Sanity checks 
    if (isempty(obj.N) || isempty(obj.M))
        error('The internal filter configuration has not been completed.');
    else  
        % Number of mixture components per plane
        J = min(obj.M * obj.N, obj.Jmax); 

        % Preallocation 
        particles = zeros(8 + 8^2, J * size(obj.planes,2));

        % Anomaly distribution initialization
        mu = linspace(0, 2*pi, J);

        % Complete particles
        for i = 1:size(obj.planes,2)
            sigma = zeros(8);
            sigma(1:end-1,1:end-1) = reshape(obj.planes(8:end,i), [7 7]);
            sigma(end,end) = deg2rad(1)^2;
            sigma = reshape(sigma, [], 1);

            % Particles as wrapped Gaussian kernels 
            particles(:,1+J*(i-1):J*i) = [repmat(obj.planes(1:7,i), 1, J); mu; repmat(sigma, 1, J)];
        end

        % Uniform generation of the weights 
        weights = repmat(obj.N/size(mu,2), 1, size(mu,2));
        weights = repmat(weights, 1, size(obj.planes,2));
    end
end