function [born_particles] = Birth(obj)
    % Constants 
    M = obj.M * obj.L + 1;

    % Uniform orbital planes birth
    born_particles(1:4,:) = obj.UniformQuat(M); 

    % Uniform orbital actions birth
    Lmin = 1; 
    Lmax = 7;
    L =  Lmin * ones(1, M) + (Lmax-Lmin) .* rand(1,M);
    emin = 0; 
    emax = 0.1; 
    e = emin * ones(1, M) + (emax-emin) .* rand(1,M);
    cos_i = repmat(-1, 1, M) + 2 .* rand(1,M);
    
    born_particles(5,:) = L; 
    born_particles(6,:) = L .* sqrt(1-e.^2); 
    born_particles(7,:) = born_particles(6,:) .* cos_i; 

    % Associated weights 
    born_particles(8,:) = repmat(1e-5/M, 1, M);
end