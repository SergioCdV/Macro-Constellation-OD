function [born_particles] = Birth(obj)
    % Constants 
    M = obj.M;

    % Preallocation 
    born_particles = zeros(8 + 8^2 + 1, size(obj.planes,2) * M);

    % Uniform orbital planes birth
    mu = linspace(0, 2*pi, M); 

    for i = 1:size(obj.planes,2)
        sigma = zeros(8);
        sigma(1:end-1,1:end-1) = reshape(obj.planes(8:end,i), [7 7]);
        sigma(end,end) = deg2rad( 10 )^2;
        sigma = reshape(sigma, [], 1);

        % Particles as wrapped Gaussian kernels 
        born_particles(1:end-1,1+M*(i-1):M*i) = [repmat(obj.planes(1:7,i), 1, M); mu; repmat(sigma, 1, M)];
    end

    % Associated weights 
    born_particles(end,:) = repmat(1/M, 1, M);
end