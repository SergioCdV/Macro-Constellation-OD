function [X] = StateEstimation(obj, samples, weights, N)
    % Configuration 
    rng(1); 

    % Extract the relevant particles 
    pruned_samples = samples(:, weights >= obj.RevThresh); 
    pruned_weights = weights(weights > obj.RevThresh);

    % Preallocation 
    X = zeros(16, N);

    % Perform K-means clustering over the quaternions and perform the state
    % estimation for the perifocal attitude
    [c, ~, index] = obj.QuatClustering(pruned_samples(1:4,:), N);

    for i = 1:N
        ID = index == i;

        % Quaternion estimation
        X(1:4,i) = obj.SteepestQuat(c(:,i), pruned_weights(1,ID), pruned_samples(1:4,ID));

        % Resampling for estimation of the action set
        [X(5:7,i), ~] = obj.Resampling(pruned_samples(5:7,ID), pruned_weights(1,ID) / sum(pruned_weights(1,ID)), 1);

        % Covariance 
        sigma = (pruned_weights(1,ID).*pruned_samples(5:7,ID)-X(5:7,i)*pruned_weights(1,ID)) * (pruned_weights(1,ID).*pruned_samples(5:7,ID)-X(5:7,i)*pruned_weights(1,ID)).';

        if (sum(ID) > 1)
            sigma = sigma/(1-sum(pruned_weights(1,ID).^2)) + 1e-6 * eye(3);
        else
            sigma = 1e-7 * eye(3);
        end
        X(8:end,i) = reshape(sigma, [], 1);
    end

    % Perform K-means clustering over the actions and perform the state
    % estimation over the plane angular velocity
%     X(1:4,:) = repmat(X(1:4,1:N), 1, N);
%     [index, a] = kmeans(pruned_samples(5:7,:).', N);
%     a = a.'; 

    % Preallocation of the covariance
%     Sigma = zeros(size(sigma,1)^2, N);
%     for i = 1:N
%         % Check the relevant particles ID 
%         ID = index == i;
% 
%         % Compute the covariance of the action set 
%         rec_samples = pruned_samples(5:7,ID);
%         a(:,i) = sum(rec_samples .* pruned_weights(ID) / sum(pruned_weights(ID), 2) ,2);
% 
%         sigma = zeros(3);
%         for j = 1:size(rec_samples,2)
%             sigma = sigma + (rec_samples(:,j)-X(5:7,i)) * (rec_samples(:,j)-X(5:7,i)).';
%         end
% 
%         if (sum(ID) > 1)
%             sigma = sigma/(sum(ID)-1) + 1e-6 * eye(3);
%         else
%             sigma = 1e-6 * eye(3);
%         end
%         X(8:end,i) = reshape(sigma, [], 1);
%     end

    % Compound the state vector 
%     for i = 1:N
%         X(5:end, 1+N*(i-1):N*i) = repmat([a(:,i); Sigma(:,i)], 1, N);
%     end
end
