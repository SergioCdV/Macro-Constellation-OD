%% Combinatorial constellation macro-determination 
% Date: 02/02/2023
% Author: Sergio Cuevas del Valle

%% Hausdorff distance %%
% This script provides the function to compute the Hausdorff distance between two sets

function [HD] = Hausdorff(P, Q, selector) 
    % Preallocation 
    D = zeros(size(P,2), size(Q,2));

    % Compute the distances
    for i = 1:size(P,2)
        for j = 1:size(Q,2)
            D(i,j) = DelaunayDistance(P(:,i), Q(:,j), selector);
        end
    end

    HD(1) = max(min(D,[],2));
    HD(2) = max(min(D.',[],2));

    HD = max(HD);
end

%% Auxiliary function
function [d] = DelaunayDistance(q2, q1, selector)
    switch (selector)
        case 1
            dq = QuaternionAlgebra.right_isoclinic(q2(1:4,1)) * QuaternionAlgebra.quaternion_inverse(q1(1:4,1)); 
            if (dq(4) < 1)
                dq = -dq;
            end
            theta = QuaternionAlgebra.MPR2Quat(1,1,dq,false);
            d = norm(theta);
        case 0
            d = norm(q2(5:7,1)-q1(5:7,1));
        case 2
            n1 = [cos(q2); sin(q2)];
            n2 = [cos(q1); sin(q1)];
            d = 1-dot(n1,n2);
    end
end