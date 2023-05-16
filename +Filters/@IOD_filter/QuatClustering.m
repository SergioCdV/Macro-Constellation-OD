

function [c, Sigma, index] = QuatClustering(obj, samples, weights)
    % Configuration 
    rng(1); 

    % Perform spherical K-means clustering over the quaternion data
%     [index, c] = kmeans(samples.', N, 'Distance', 'cosine');
%     c = (c./sqrt(dot(c,c,2))).';
%     Sigma = [];

    % Perform density-based clustering 
    epsilon = 0.01;
    minpts = max(1, round( 0.1 * size(samples,2) ));

    D = pdist(samples, "cosine");

    if (size(D,1))
        index = dbscan(D, epsilon, minpts, "Distance", "precomputed");
    
        c = zeros(size(samples,1), max(index));
        for i = 1:max(index)
            sel = samples(1:4, index == i);
            c(:,i) = sel(:,1);
        end
        
        if (any(index < 0))
            c = [c samples(:, index == -1)];
        end
    else
        c = [];
        index = [];
    end

    Sigma = [];
end