

function [c, Sigma, index] = QuatClustering(obj, weights, samples)
    % Density between quaternions 
%     mu = mean(weights);
%     sigma = std(weights);
%     samples = samples(:, abs(weights-mu) < sigma); 
    d = real( pdist(samples.', @personal_distance) );
    D = squareform( d );

    % Perform density-based clustering 
    epsilon = deg2rad( 15 );
    minpts = size(samples,1) + 1;

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

%% Auxiliary functions 
% NANEUCDIST Euclidean distance ignoring coordinates with NaNs
function [D] = personal_distance(XI,XJ)  
    % Correction for missing coordinates
    d = XI * XJ.';
    D = acos(2*d.^2-1);
end

%% Auxiliary code 

    % Perform spherical K-means clustering over the quaternion data
%     [index, c] = kmeans(samples.', N, 'Distance', 'cosine');
%     c = (c./sqrt(dot(c,c,2))).';
%     Sigma = [];