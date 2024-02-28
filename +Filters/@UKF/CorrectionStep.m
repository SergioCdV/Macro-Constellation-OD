
function [X, P, Pmeas, Y] = CorrectionStep(obj, sigma, X, P, z)
    % Measurement prediction 
    y = feval(obj.ObservationModel, sigma);
    Y = measurements_prediction(obj, y);
    
    % State and covariance prediction 
    switch (obj.Algorithm)
        case 'UKF-A'
            [X, P, Pmeas] = UKFA_correction(obj, sigma, X, P, y, Y, z);
            
        case 'UKF-S'
            [X, P, Sy] = UKFS_correction(obj, sigma, X, P, y, Y, z);
            Pmeas = Sy*Sy.';
    end
end