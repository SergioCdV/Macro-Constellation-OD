
function [State, Sigma, Pmeas, y] = CorrectionStep(obj, sigma, State, Sigma, z)
    % Compute a mean quaternion state
    QuatState = State(end-3:end,:);
    X = [sigma(end-3:end,:); sigma(4:end-4,:)];

    % Measurement prediction 
    Y = feval(obj.ObservationModel, X);
    y = measurements_prediction(obj, Y);
    
    % State and covariance prediction 
    switch (obj.Algorithm)
        case 'UKF-A'
            [x, Sigma, Pmeas] = UKFA_correction(obj, sigma(1:end-4,:), State(1:end-4,:), Sigma, Y, y, z);
            
        case 'UKF-S'
            [x, Sigma, Sy] = UKFS_correction(obj, sigma(1:end-4,:), State(1:end-4,:), Sigma, Y, y, z);
            Pmeas = Sy*Sy.';
    end

    dq = QuaternionAlgebra.MPR2Quat(obj.a, obj.f, x(1:3,:), true);
    State = [x; QuaternionAlgebra.right_isoclinic(dq) * QuatState];
end