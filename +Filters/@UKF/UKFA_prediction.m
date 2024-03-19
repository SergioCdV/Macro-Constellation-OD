%% Constellation macro-orbit determination %%
% Date: 19/03/2024

%% Propagation step in additive UKF
% This script provides the implementation of the propagation step in an additive UKF

function [X, P] = UKFA_prediction(obj, sigma)
    % State prediction
    X = sum(obj.W(1,:) .* sigma, 2);

    % Covariance prediction
    aux_state = sigma - X;
    P = obj.Q + ( obj.W(2,:) .* aux_state) * aux_state.';
end