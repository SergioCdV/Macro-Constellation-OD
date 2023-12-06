
% Jose-Luis Blanco (2023). Resampling methods for particle filtering (https://www.mathworks.com/matlabcentral/fileexchange/24968-resampling-methods-for-particle-filtering), MATLAB Central File Exchange. Recuperado May 13, 2023.

function [particles, weights] = Resampling(obj, particles, weights, N)
    switch (obj.ResamplingMethod)
        case 'Residual'
            [particles, weights] = ResamplingResidual(particles, weights, N);
        case 'Multinomial'
            [particles, weights] = ResamplingMultinomial(particles, weights, N);
        case 'Systematic'
            [particles, weights] = ResamplingSystematic(particles, weights, N);
        case 'Stratified'
            [particles, weights] = ResamplingStratified(particles, weights, N);
        otherwise
            error('No valid resampling method has been selected.');
    end
end

%% Auxiliary functions 
% Residual resampling 
function [particles, weights] = ResamplingResidual(particles, weights, N)
    M = N;

    % Deterministic part
    Ns = floor(M .* weights);                           % "Repetition counts" (plus the random part, later on):
    weights = weights - Ns/N;
    old_particles = particles;
    particles = zeros(size(old_particles,1), M);

    n = 1; 
    for i = 1:length(weights)
        for h = 1:weights(i)
            n = n+1;
            particles(:,n) = old_particles(:,i);
        end
    end
    Ns = n;

    % Complete the set using multinomial resampling
    weights = weights * N / (N-Ns);
    [particles(:,Ns+1:N), ~] = ResamplingMultinomial(particles, weights, N-Ns);

    % Sample the particles 
    weights = repmat(1/size(particles,2),1,size(particles,2));
end

% Multinomial resampling 
function [particles, weights] = ResamplingMultinomial(particles, weights, N)
    % Preallocation 
    old_particles = particles;
    particles = zeros(size(particles,1), N);
    Q = cumsum(weights);

    % Get the indices
    i = 1;
    while (i <= N)
        sampl = rand();
        j = 1;
        while (Q(j) < sampl)
            j = j + 1;
        end
        particles(:,i) = old_particles(:,j);
        i = i + 1;
    end

    % Sample the particles 
    weights = repmat(1/size(particles,2),1,size(particles,2));
end

% Systematic resampling 
function [particles, weights] = ResamplingSystematic(particles, weights, N)
    % Preallocation 
    old_particles = particles;
    particles = zeros(size(particles,1), N);
    Q = cumsum(weights);

    T = linspace(0,1-1/N,N) + rand()/N;
    T(N+1) = 1;

    % Get the indices
    i = 1;
    j = 1;
    while (i <= N && j)
        if (T(i) >= Q(j))
            j = j + 1;

            if (j >= size(Q,2))
                j = size(Q,2);
            end
        else
            particles(:,i) = old_particles(:,j);
            i = i + 1;
        end
    end

    % Sample the particles 
    weights = repmat(1/size(particles,2),1,size(particles,2));
end

% Stratified resampling 
function [particles, weights] = ResamplingStratified(particles, weights, N)
    % Preallocation 
    old_particles = particles;
    particles = zeros(size(particles,1), N);
    Q = cumsum(weights);

    T = zeros(1,N);
    for i = 1:N
        T(i) = rand(1,1)/N + (i-1)/N;
    end
    T(N+1) = 1;

    % Get the indices
    i = 1;
    j = 1;
    while (i <= N)
        if (T(i) >= Q(j))
            j = j + 1;
        else
            particles(:,i) = old_particles(:,j);
            i = i + 1;
        end
    end

    % Sample the particles 
    weights = repmat(1/size(particles,2),1,size(particles,2));
end
