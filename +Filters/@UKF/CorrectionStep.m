
function [State, Sigma, Pmeas, y] = CorrectionStep(obj, sigma, State, Sigma, z)
    % Measurement prediction 
    y = feval(obj.ObservationModel, sigma);
    Y = measurements_prediction(obj, y);
    
    % State and covariance prediction 
    switch (obj.Algorithm)
        case 'UKF-A'
            [State, Sigma, Pmeas] = UKFA_correction(obj, sigma, State, Sigma, y, Y, z);
            
        case 'UKF-S'
            [State, Sigma, Sy] = UKFS_correction(obj, sigma, State, Sigma, y, Y, z);
            Pmeas = Sy*Sy.';
    end
end