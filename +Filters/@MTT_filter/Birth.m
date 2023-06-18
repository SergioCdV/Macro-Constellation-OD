function [born_particles] = Birth(obj)
    % Constants 
    M = obj.M;
    pos = 8;

    % Preallocation 
    born_particles = zeros(pos + (pos-1)^2 + 1, size(obj.planes,2) * M);

    % Uniform orbital planes birth
    mu = linspace(0, 2*pi, M); 

    for i = 1:size(obj.planes,2)
        sigma = zeros(pos-1);
        sigma(1:end-1,1:end-1) = reshape(obj.planes(pos:end,i), [pos-2 pos-2]);
        sigma(end,end) = deg2rad( 10 );
        sigma = reshape(sigma, [], 1);

        % Particles as wrapped Gaussian kernels 
        born_particles(1:end-1,1+M*(i-1):M*i) = [repmat(obj.planes(1:pos-1,i), 1, M); mu; repmat(sigma, 1, M)];
    end

    % Associated weights 
    born_particles(end,:) = repmat(1/M, 1, M);
end