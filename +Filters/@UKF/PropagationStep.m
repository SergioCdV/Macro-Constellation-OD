%% Constellation macro-orbit determination %%
% Date: 19/03/2024

%% Propagation step UKF
% This script provides the implementation of the generic propagation step in a UKF

function [sigma, X, P] = PropagationStep(obj, time_step)
    % Generate sigma points 
    sigma = sigma_points(obj, obj.State, obj.Sigma);

    % Propagation of sigma points 
    sigma = feval(obj.StateModel, sigma, time_step);

    % State and covariance prediction 
    switch (obj.Algorithm)
        case 'UKF-A'
            [X, P] = UKFA_prediction(obj, sigma);
        case 'UKF-S'
            [X, P] = UKFS_prediction(obj, sigma);
    end
end