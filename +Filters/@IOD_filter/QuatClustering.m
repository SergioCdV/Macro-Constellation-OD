

function [c, Sigma, index] = QuatClustering(obj, samples, N)
    % Configuration 
    rng(1); 

    % Perform K-means clustering over the quaternion data
    [index, c] = kmeans(samples.', N, 'Distance', 'cosine');
    Sigma = [];
    c = c.';
end