

function [sigma, State, Sigma] = PropagationStep(obj, time_step)
    % Generate sigma points 
    sigma = sigma_points(obj, obj.State, obj.Sigma);

    % Propagation of sigma points 
    sigma = feval(obj.StateModel, sigma, time_step);

    % State and covariance prediction 
    switch (obj.Algorithm)
        case 'UKF-A'
            [State, Sigma] = UKFA_prediction(obj, sigma);
        case 'UKF-S'
            [State, Sigma] = UKFS_prediction(obj, sigma);
    end
end