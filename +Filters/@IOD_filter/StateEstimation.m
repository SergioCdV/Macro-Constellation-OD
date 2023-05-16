function [X] = StateEstimation(obj, samples, weights, T)
    % Configuration 
    rng(1); 

    % Extract the relevant particles 
    pruned_samples = samples(:, weights >= obj.RevThresh); 
    pruned_weights = weights(weights > obj.RevThresh);

    % Perform K-means clustering over the quaternions and perform the state estimation for the perifocal attitude
    [c, ~, index] = obj.QuatClustering(pruned_samples(1:4,:));

    % Estimate the number of planes based on the quaternion distribution
    N = size(c,2);

    % Preallocation 
    X = zeros(16, N);
    pos = 0;

    for i = 1:max(index)
        ID = index == i;
        pos = i;

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

    if (any(index < 0))
        ID = index == -1;
        
        for i = pos+1:max(index)+sum(ID)
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
    end


end
