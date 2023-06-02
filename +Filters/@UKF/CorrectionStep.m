
function [State, Sigma, Pmeas, y] = CorrectionStep(obj, sigma, State, Sigma, z)
    % Measurement prediction 
    Y = feval(obj.ObservationModel, sigma);
    y = measurements_prediction(obj, Y);
    
    % State and covariance prediction 
    switch (obj.Algorithm)
        case 'UKF-A'
            [State, Sigma, Pmeas] = UKFA_correction(obj, sigma, State, Sigma, Y, y, z);
            
        case 'UKF-S'
            [State, Sigma, Sy] = UKFS_correction(obj, sigma, State, Sigma, Y, y, z);
            Pmeas = Sy*Sy.';
    end
end