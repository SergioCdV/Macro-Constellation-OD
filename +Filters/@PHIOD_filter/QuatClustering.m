

function [c, Sigma, index] = QuatClustering(obj, weights, samples)
    % Density between quaternions 
    % d = real( pdist(samples.', @personal_distance) );
    samples(1:4,:) = samples(1:4,:) ./ sqrt(dot(samples(1:4,:), samples(1:4,:), 1));
    D = real( personal_distance(samples, samples)) ; %squareform( d );
 
    % Perform density-based clustering 
    epsilon = deg2rad( 15 );
    minpts = 3;%size(samples,1) + 1;

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

    [c, ia, ~] = unique(c.', 'rows');
    c = c.';
    index = index(ia);

    Sigma = [];
end

%% Auxiliary functions 
function [D] = personal_distance(XI,XJ)  
    % Correction for missing coordinates
    D = zeros(size(XI,2), size(XJ,2));
    for i = 1:size(XI, 2)
        inv = QuaternionAlgebra.quaternion_inverse( XI(1:4,i) ); 

        for j = 1:size(XJ,2)
            aux = QuaternionAlgebra.right_isoclinic(inv) * XJ(1:4,j);
            D(i,j) = 2 * acos(aux(4));
        end
    end
end