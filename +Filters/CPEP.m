%% Combinatorial constellation macro-determination 
% Date: 02/02/2023
% Author: Sergio Cuevas del Valle

%% Circular Error Probable %%
% This script provides the function to compute the CPEP for a given sequence of multi-target estimates

function [CPEP] = CPEP(N, N_hat, X, X_hat) 
    % Preallocation
    rho = zeros(1,N);
    outlier = 0;
    counter = 0;
    
    for j = 1:N
        for k = 1:N_hat
            d = DelaunayDistance(X(:,j), X_hat(:,k));
            if (d > r)
               outlier = outlier+1; 
            end
            counter = counter + 1; 
        end
        rho(j) = outlier/counter;
    end

    CPEP = sum(rho) / N;
end