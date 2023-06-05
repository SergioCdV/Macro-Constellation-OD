
function [X, P] = UKFA_prediction(obj, sigma)
    % State prediction
    X = sum(obj.W(1,:).*sigma,2);

    % Covariance prediction
    P = obj.Q;
    for i = 1:size(sigma,2)
        aux = sigma(:,i)-X;
        P = P + obj.W(2,i) * (aux * aux.');
    end
end