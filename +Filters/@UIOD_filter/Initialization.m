

function [particles, weights] = Initialization(obj)
    % Sanity checks 
    if (isempty(obj.N) || isempty(obj.M))
        error('The internal filter configuration has not been completed.');
    else  
        % Number of mixture components per plane
        J = max(obj.Jmax, obj.M); 
        pos = 7;

        % Preallocation 
        particles = zeros(pos + (pos-1)^2, J);

        % Quaternion distribution initialization
        particles(1:4,:) = QuaternionAlgebra.UniformQuat(J);

        % Action set 
        Lmin = 1; 
        Lmax = 1.5;
        L =  Lmin * ones(1,J) + (Lmax-Lmin) .* rand(1,J);
        emin = 0; 
        emax = 0.01; 
        e = emin * ones(1, J) + (emax-emin) .* rand(1,J);
        cos_i = repmat(-1, 1, J) + 2 .* rand(1,J);
        
        particles(5,:) = L; 
        particles(6,:) = L .* sqrt(1-e.^2); 
        particles(7,:) = particles(6,:) .* cos_i; 

%         mu = [1.04; 1e-3; cos(deg2rad(45))];
%         particles(5:7,:) = repmat(mu, 1, J);
%         particles(6,:) = particles(5,:) .* sqrt(1-particles(6,:).^2); 
%         particles(7,:) = particles(6,:) .* particles(7,:);
%         particles(1:4,:) = repmat([0.2164 0 0 0.9763].', 1, size(particles,2));

        % Covariance matrix
        sigma = 1e-3 * eye(pos-1);
        sigma = reshape(sigma, [], 1);

        % Particles as wrapped Gaussian kernels 
        particles(pos+1:end,:) = repmat(sigma, 1, J);

        % Uniform generation of the weights 
        weights = repmat(obj.N/size(particles,2), 1, size(particles,2));
    end
end