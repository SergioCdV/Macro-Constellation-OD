function [born_particles] = Birth(obj)
    % Constants 
    J = obj.M;          % Number of birth particles
    pos = 7;            % State dimensions

    % Preallocation 
    born_particles = zeros(pos + (pos-1)^2 + 1, J);

    % Uniform orbital planes birth
    born_particles(1:4,:) = QuaternionAlgebra.UniformQuat(J);

    % Action set 
    Lmin = 1; 
    Lmax = 7;
    L =  Lmin * ones(1,J) + (Lmax-Lmin) .* rand(1,J);
    emin = 0; 
    emax = 0.1; 
    e = emin * ones(1, J) + (emax-emin) .* rand(1,J);
    cos_i = repmat(-1, 1, J) + 2 .* rand(1,J);
    
    born_particles(5,:) = L; 
    born_particles(6,:) = L .* sqrt(1-e.^2); 
    born_particles(7,:) = born_particles(6,:) .* cos_i; 

    % Covariances
    sigma = 1e-7 * eye(pos-1);
    sigma = reshape(sigma, [], 1);

    % Particles as wrapped Gaussian kernels 
    born_particles(pos+1:end-1,:) = repmat(sigma, 1, J);

    % Associated weights 
    born_particles(end,:) = repmat(1/J, 1, J);
end