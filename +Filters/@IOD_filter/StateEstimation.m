function [X] = StateEstimation(obj, samples, weights, N)
    % Configuration 
    rng(1); 

    % Perform K-means clustering over the quaternions
    [c, ~, index] = obj.QuatClustering(samples(1:4,:), N);

    % Preallocation 
    X = [c; zeros(12,size(c,2))];
        
    % State estimation
    for i = 1:size(c,2)
        % Assemble the appropriate particles per cluster 
        ID = index(index == i);

        % Compute the new quaternions
        X(1:4,i) = SteepestQuat(X(1:4,i), weights(1,ID), samples(1:4,ID));

        % Compute the associate action set
        X(5:7,i) = sum( weights(1,ID) .* samples(5:7, ID), 2);

        % Compute the covariance of the action set 
        sigma = (samples(5:7, ID)-X(5:7,i)) * (samples(5:7, ID)-X(5:7,i)).';
        X(8:end,i) = reshape(sigma/length(ID), [], 1);
    end
end

%% Auxiliary function 
function [q] = SteepestQuat(q0, w, dq)
    % Set up the Newton Rhapson method 
    maxIter = 100; 
    iter = 1; 
    tol = 1e-3; 
    GoOn = true; 
    q = q0;

    while (GoOn && iter < maxIter)
        % Compute the update term 
        dQ = zeros(4,1);
        for i = 1:size(dq,2)
            dQ = dQ + w(i) * QuaternionAlgebra.log_map(dq(:,i),q);
        end

        qn = QuaternionAlgebra.exp_map(dQ,q);
        res = q.' * qn;
        q = qn;

        % Convergence analysis
        if (abs(1-res) < tol)
            GoOn = false;
        else
            iter = iter + 1;
        end
    end
end