
function [X, P] = UKFA_prediction(obj, sigma)
    % State prediction
    X = sum(obj.W(1,:).*sigma,2);

    % Covariance prediction
    P = (sigma-X)*diag(obj.W(2,:))*(sigma-X).'+obj.Q;
end