function [X] = StateEstimation(obj, samples, weights, N)
    % Configuration 
    rng(1); 

    % Perform K-means clustering over the quaternions
    [c, ~, index] = obj.QuatClustering(samples(1:4,:), N);

    % Preallocation 
    X = [c; zeros(12,size(c,2))];

    % State estimation
    for i = 1:size(X,2)
        % Assemble the appropriate particles per cluster 
        ID = index == i;

        % Compute the new quaternions
        X(1:4,i) = obj.SteepestQuat(X(1:4,i), weights(1,ID), samples(1:4,ID));

        % Compute the associate action set
        weights(1,ID) = weights(1,ID) / sum( weights(1,ID) );
        X(5:7,i) = sum( weights(1,ID) .* samples(5:7, ID), 2);

        % Compute the covariance of the action set 
        sigma = zeros(3);
        rec_samples = samples(5:7,ID);

        for j = 1:size(rec_samples,2)
            sigma = sigma + (rec_samples(:,j)-X(5:7,i)) * (rec_samples(:,j)-X(5:7,i)).';
        end
        sigma = sigma/(sum(ID)-1) + 1e-6 * eye(3);
        X(8:end,i) = reshape(sigma, [], 1);
    end
end
