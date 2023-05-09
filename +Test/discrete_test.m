%% Combinatorial constellation macro-determination 
% Date: 27/10/2022
% Author: Sergio Cuevas del Valle

%% Input data 
N = 3;                          % Number of spacecraft 
dtheta = deg2rad(17);           % Slot resolution

M = floor(2*pi/dtheta);         % Number of slots
slots = linspace(0,2*pi,M);     % Anomaly slots 

%% PFM computation
% Compute the ring anomaly distributions and the associated constellation discrete PMF
[alpha] = possible_distributions(dtheta, N);
[X] = discrete_pmf(dtheta, N);

%% Auxiliary functions 
% Discrete probability density function 
function [X] = discrete_pmf(dtheta, N)
    % Compute the constellations 
    Constellations = possible_distributions(dtheta, N);

    % Apply normalization 
    Constellations = Constellations / N; 
    X = Constellations / size(Constellations,1);
end

% Compute the possible distribution sets for each orbital plane
function [alpha, M] = possible_distributions(dtheta, N)
    % Compute the possible set of constellations 
    M = floor(2*pi/dtheta);

    % Sanity check on the slots 
    if (M < N)
        warning('Number of slots is not sufficient'); 
        M = N; 
    end

    % Compute the associated Gray code 
    alpha = gray(M);

    % Restrict the code by particle conservation 
    alpha = alpha(sum(alpha,2) == N,:);
end

% Compute the Grat code of dimensions M recursively 
function [alpha] = gray(M)
    % Branch 
    if (M == 1)
        alpha = [0 1].';
    else
        beta = gray(M-1);
        gamma = flip(beta,1);
        beta = reshape(beta.', 1, size(beta,1)*size(beta,2));
        gamma = reshape(gamma.', 1, size(gamma,1)*size(gamma,2));

        % Go for the recursion
        alpha = zeros(1,2^M*M);
        
        for i = 1:2^(M-1)
            alpha(1+(M)*(i-1)) = 0;
            alpha(2^M*M/2+1+(M)*(i-1)) = 1;
        end

        for i = 1:2^(M-1)
            alpha(2+(M)*(i-1):2+(M)*(i-1)+(M-1)-1) = beta(1+(M-1)*(i-1):(M-1)*i);
            alpha(2^M*M/2+2+(M)*(i-1):2^M*M/2+2+(M)*(i-1)+(M-1)-1) = gamma(1+(M-1)*(i-1):(M-1)*i);
        end

        alpha = reshape(alpha, [M, 2^M]).';
    end
end
