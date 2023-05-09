%% Combinatorial constellation macro-determination 
% Date: 02/02/2023
% Author: Sergio Cuevas del Valle

%% Hausdorff distance %%
% This script provides the function to compute the Hausdorff distance between two sets

function [HD] = Hausdorff(P, Q) 
    % Preallocation 
    D = zeros(size(P,1), size(Q,1));

    % Compute the distances
    for i = 1:size(P,1)
        for j = 1:size(Q,1)
            D(i,j) = DelaunayDistance(P(:,i), Q(:,j));
        end
    end

    HD(1) = max(min(D,[],2));
    HD(2) = max(min(D.',[],2));

    HD = max(HD);
end