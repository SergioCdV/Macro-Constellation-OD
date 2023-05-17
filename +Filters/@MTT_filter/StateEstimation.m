function [X] = StateEstimation(obj, samples, weights, N)
    % Configuration 
    rng(1); 

    % Extract the relevant particles through resampling
    T = sum(weights);
    [pruned_samples, ~] = obj.Resampling(samples, weights / T, size(samples,2)); 
    pruned_weights = weights;

    % Perform K-means clustering over the anomaly
    N = min(N, size(pruned_samples,2));

    if (N)
        [index, X(:,1)] = kmeans(mod(pruned_samples(1,:), 2*pi), N);
        X = X.'; 
    
        for i = 1:max(index)
            ID = index == i;
    
            pruned_weights(1,ID) = pruned_weights(1,ID) / sum(pruned_weights(1,ID));
    
            % Covariance 
            sigma = (pruned_weights(1,ID).*pruned_samples(1,ID)-X(1,i)*pruned_weights(1,ID)) * (pruned_weights(1,ID).*pruned_samples(1,ID)-X(1,i)*pruned_weights(1,ID)).';
    
            if (sum(ID) > 1)
                sigma = sigma/(1-sum(pruned_weights(1,ID).^2)) + 1e-7;
            else
                sigma = 1e-7;
            end
            X(2:end,i) = reshape(sigma, [], 1);
        end
    else
        X = [];
    end

end
