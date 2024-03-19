%% Constellation macro-orbit determination %%
% Date: 19/03/2024

%% Correction step in square root UKF
% This script provides the implementation of the correction step in an square root additive UKF

function [X, S] = UKFS_prediction(obj, sigma)
    % State prediction
    X = sum(obj.W(1,:) .* sigma, 2);

    % Covariance prediction
    res = sqrt( abs(obj.W(2,:)) ) .* ( sigma - X );
    [~, S] = qr( [res(:,2:end) obj.Q].', 0 );

    if (obj.W(2,1) < 0)
        S = cholupdate(S, res(:,1), '-');
    else
        S = cholupdate(S, res(:,1), '+');
    end
end