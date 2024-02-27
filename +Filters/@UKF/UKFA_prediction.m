
function [X, P] = UKFA_prediction(obj, sigma)
    % State prediction
    X = sum(obj.W(1,:) .* sigma, 2);

    % Covariance prediction
    aux_state = sigma - X;
    P = obj.Q + ( obj.W(2,:) .* aux_state) * aux_state.';
end