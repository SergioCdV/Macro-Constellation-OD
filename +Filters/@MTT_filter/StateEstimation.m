function [X] = StateEstimation(obj, samples, weights, N)
    % Configuration 
    rng(1); 
    pos = 8;

    % Extract the relevant particles through resampling
    [pruned_samples, ~] = obj.Resampling(samples, weights / sum(weights), size(samples,2)); 
    pruned_weights = weights;

    % Perform K-means clustering over the anomaly
    N = min(N, size(pruned_samples,2));

    if (N)
        X = zeros(pos + pos^2, N);
        [index, c] = kmeans(mod(pruned_samples(pos,:).', 2*pi), N);
        X(pos,:) = c.'; 
    
        for i = 1:max(index)
            ID = index == i;
    
            pruned_weights(1,ID) = pruned_weights(1,ID) / sum(pruned_weights(1,ID));
    
            % Covariance 
            sigma = (pruned_weights(1,ID).*pruned_samples(pos,ID)-X(1,i)*pruned_weights(1,ID)) * (pruned_weights(1,ID).*pruned_samples(pos,ID)-X(1,i)*pruned_weights(1,ID)).';
    
            if (sum(ID) > 1)
                sigma = sigma/(1-sum(pruned_weights(1,ID).^2)) + 1e-7;
            else
                sigma = 1e-7;
            end
            X(end,i) = reshape(sigma, [], 1);
        end
        X([1:pos pos+1:end],:) = repmat(pruned_samples([1:pos pos+1:end],1), 1, N);
    else
        X = [];
    end

end
