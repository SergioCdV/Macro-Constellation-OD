function [X] = StateEstimation(obj, samples, weights, N)
    % Configuration  
    pos = 8;

    % Extract the relevant particles through resampling
    [pruned_samples, ~] = obj.Resampling(samples, weights / sum(weights), size(samples,2)); 
    pruned_weights = weights;

    % Perform K-means clustering over the anomaly
    N = min(N, size(pruned_samples,2));

    if (N)
        X = zeros(pos + (pos-1)^2, N);
        v = [cos(pruned_samples(pos,:).') sin(pruned_samples(pos,:).')];
        [index, c] = kmeans(v, N, 'Distance','cosine');
        X(pos,:) = atan2(c(:,2), c(:,1)).'; 

        for i = 1:max(index)
            ID = index == i;
    
            pruned_weights(1,ID) = pruned_weights(1,ID) / sum(pruned_weights(1,ID));
    
            % Covariance 
            sigma = (pruned_weights(1,ID).*pruned_samples(pos,ID)-X(1,i)*pruned_weights(1,ID)) * (pruned_weights(1,ID).*pruned_samples(pos,ID)-X(1,i)*pruned_weights(1,ID)).';
    
            if (sum(ID) > 1)
                sigma = sigma/(1-sum(pruned_weights(1,ID).^2)) + obj.PD_tol;
            else
                sigma = obj.PD_tol;
            end
            X(end,i) = reshape(sigma, [], 1);
        end
        
        X([1:pos-1 pos+1:end-1],:) = repmat(pruned_samples([1:pos-1 pos+1:end-1],1), 1, N);

    else
        X = [];
    end

end
