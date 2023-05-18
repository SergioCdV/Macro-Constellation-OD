
function [State, Sigma, Pyy] = UKFA_correction(obj, sigma, State, Sigma, Y, y, z)
    % Covariances matrices
    Pyy = (Y-y)*diag(obj.W(2,:))*(Y-y).'+obj.R;        
    Pxy = (sigma-State)*diag(obj.W(2,:))*(Y-y).';
    
    % Kalman gain
    K = Pxy*(Pyy^(-1));

    % Update
    State = State+K*(z-y);
    Sigma = Sigma-K*Pyy*K.';
end