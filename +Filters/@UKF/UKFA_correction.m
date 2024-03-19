%% Constellation macro-orbit determination %%
% Date: 19/03/2024

%% Correction step in additive UKF
% This script provides the implementation of the correction step in an additive UKF

function [State, Sigma, Pyy] = UKFA_correction(obj, sigma, State, Sigma, y, Y, z)
    % Update
    if (~any(isnan(y)))
        % Covariances matrices        
        aux_state = sigma - State;
        aux_meas = y - Y;
        Pxy = (obj.W(2,:) .* aux_state) * aux_meas.';             % Cross correlation matrix
        Pyy = obj.R + (obj.W(2,:) .* aux_meas) * aux_meas.';      % Measurement covariance
        
        % Kalman gain
        K = Pxy * ( Pyy \ eye(size(Pyy)) );

        State = State+K*(z-Y);     
        Sigma = Sigma-Pxy*K.'-K*Pxy.'+K*Pyy*K.';    % Joseph update 
    else
        Pyy = [];
    end
end