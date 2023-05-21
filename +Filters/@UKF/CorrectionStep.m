
function [State, Sigma, Pmeas, Y] = CorrectionStep(obj, sigma, State, Sigma, z)
    % Measurement prediction 
    y = feval(obj.ObservationModel, sigma);
    Y = measurements_prediction(obj, y);
    
    % State and covariance prediction 
    switch (obj.Algorithm)
        case 'UKF-A'
            [State, Sigma, Pmeas] = UKFA_correction(obj, sigma, State, Sigma, y, Y, z);
            Sigma = 0.5 * (Sigma + Sigma.') + 1E-6 * eye(size(Sigma,1));
            
        case 'UKF-S'
            [State, Sigma, Sy] = UKFS_correction(obj, sigma, State, Sigma, y, Y, z);
            Pmeas = Sy*Sy.';
    end
end