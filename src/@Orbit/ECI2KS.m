%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 17/11/22

%% State mapping %%
% Function to compute the mapping between the Cartesian state vector and
% the u space

% Inputs: - array x, either the Cartesian of the u-space state vector
%         - boolean direction, to determine the sense of the transformation

% Outputs: - array S, the transformed state vector

function [S] = ECI2KS(x, direction)
    % Compute the mapping
    if (direction)
        % Preallocation 
        S = zeros(size(x,1),8); 

        % Transformation from the Cartesian space to the U space
        for i = 1:size(x,1)
            S(i,1:4) = u_mapping(x(i,1:3));              % Position space transformation
            L = KS_matrix(S(i,1:4));                     % KS matrix
            S(i,5:8) = (1/2)*L.'*[x(i,4:6); 0];          % Velocity space transformation
        end
    else
        % Preallocation 
        S = zeros(size(x,1),6);

        % Transformation from the U space to the Cartesian space
        for i = 1:size(x,1)
            L = KS_matrix(x(i,1:4));                                 % KS matrix
            aux = L*x(i,1:4);                                        % Position space transformation
            S(i,1:3) = aux(1:3);                                     % Position space transformation
            aux = 2/dot(S(i,1:4),S(i,1:4))*L*x(i,6:9);               % Velocity space transformation
            S(i,4:6) = aux(1:3);                                     % Velocity space transformation
        end
    end
end

%% Auxiliary functions 
% Posotion to U space mapping
% Inputs: - array r, the Cartesian position vector, an m x r array
% Outputs: - array r, the KS position vector, an m x r array
function [u] = u_mapping(r)
    % Pre-allocation 
    u = zeros(size(r,1),4); 

    % Compute the mapping to the Hopf fibre 
    for i = 1:size(r,1)
        if (r(i,1) >= 0)
            u(i,4) = 0;
            u(i,1) = ((r(i,1)+norm(r(i,:)))/2)^(1/2);
            u(i,2) = (1/2)*(r(i,2)/u(i,1));
            u(i,3) = (1/2)*(r(i,3)/u(i,1));
        elsei
            u(i,3) = 0;
            u(i,2) = ((norm(r(i,:))-r(i,1))/2)^(1/2);
            u(i,1) = (1/2)*(r(i,2)/u(i,2));
            u(i,4) = (1/2)*(r(i,3)/u(i,2));
        end
    end
end