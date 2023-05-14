function [X] = StateEstimation(obj, samples, weights, N)
    % Configuration 
    rng(1); 

    % Extract the relevant particles 
    pruned_samples = samples(:, weights >= obj.RevThresh); 
    pruned_weights = weights(weights > obj.RevThresh);

    % Preallocation 
    X = zeros(16, N^2);

    % Perform K-means clustering over the quaternions and perform the state
    % estimation for the perifocal attitude
    [c, ~, index] = obj.QuatClustering(pruned_samples(1:4,:), N);
    for i = 1:N
        ID = index == i;
        X(1:4,i) = obj.SteepestQuat(c(:,i), pruned_weights(1,ID), pruned_samples(1:4,ID));
    end

    X(1:4,:) = repmat(X(1:4,1:N), 1, N);

    % Perform K-means clustering over the actions and perform the state
    % estimation over the plane angular velocity
    [index, a] = kmeans(pruned_samples(5:7,:).', N);
    a = a.'; 

    % Preallocation of the covariance
    sigma = zeros(size(a,1));
    Sigma = zeros(size(a,1)^2, N);
    for i = 1:N
        % Check the relevant particles ID 
        ID = index == i;

        % Compute the covariance of the action set 
        rec_samples = pruned_samples(5:7,ID);
        a(:,i) = sum(rec_samples .* pruned_weights(ID) / sum(pruned_weights(ID), 2) ,2);

        for j = 1:size(rec_samples,2)
            sigma = sigma + (rec_samples(:,j)-a(:,i)) * (rec_samples(:,j)-a(:,i)).';
        end
        sigma = sigma/(sum(ID)-1) + 1e-6 * eye(3);
        Sigma(:,i) = reshape(sigma, [], 1);
    end

    % Compound the state vector 
    for i = 1:N
        X(5:end, 1+N*(i-1):N*i) = repmat([a(:,i); Sigma(:,i)], 1, N);
    end
end
