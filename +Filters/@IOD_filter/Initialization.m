

function [particles, weights] = Initialization(obj)
    % Sanity checks 
    if (isempty(obj.N) || isempty(obj.L) || isempty(obj.M))
        error('The internal filter configuration has not been completed.');
    else      
        % Initial uniform distribution on the quaternion sphere 
        M = obj.N * (obj.L * obj.M + 1);
        particles(1:4,:) = obj.UniformQuat(M);

        % Generate the initial uniform distribution on the action space
        Lmin = 1; 
        Lmax = 7;
        L =  Lmin * ones(1, M) + (Lmax-Lmin) .* rand(1,M);
        emin = 0; 
        emax = 0.1; 
        e = emin * ones(1, M) + (emax-emin) .* rand(1,M);
        cos_i = repmat(-1, 1, M) + 2 .* rand(1,M);
        
        particles(5,:) = L; 
        particles(6,:) = L .* sqrt(1-e.^2); 
        particles(7,:) = particles(6,:) .* cos_i; 

        mu = [1.08; 0.001; cos(deg2rad(50))];
        particles(5:7,:) = repmat(mu, 1, M);

        % Uniform generation of the weights 
        weights = repmat(1/size(particles,2), 1, size(particles,2));
    end
end