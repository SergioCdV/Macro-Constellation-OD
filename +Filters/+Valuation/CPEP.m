%% Combinatorial constellation macro-determination 
% Date: 02/02/2023
% Author: Sergio Cuevas del Valle

%% Circular Error Probable %%
% This script provides the function to compute the CPEP for a given sequence of multi-target estimates

function [CPEP] = CPEP(N, N_hat, X_hat, X, r, selector) 
    % Preallocation
    rho = zeros(1,N);
    outlier = 0;
    counter = 0;
    
    for j = 1:N
        for k = 1:N_hat
            d = DelaunayDistance(X(:,j), X_hat(:,k), selector);
            if (d > r)
               outlier = outlier+1; 
            end
            counter = counter + 1; 
        end
        rho(j) = outlier/counter;
    end

    CPEP = sum(rho) / N;
end

%% Auxiliary function
function [d] = DelaunayDistance(q2, q1, selector)
    switch (selector)
        case 1
            dq = QuaternionAlgebra.right_isoclinic(q2(1:4,1)) * QuaternionAlgebra.quaternion_inverse(q1(1:4,1)); 
            theta = QuaternionAlgebra.MPR2Quat(1,1,dq,false);
            d = norm(theta);
        case 0
            d = norm(q2(5:7,1)-q1(5:7,1));
        case 2
            d = norm([cos(q2); sin(q2)] - [cos(q1); sin(q1)]);
    end
end