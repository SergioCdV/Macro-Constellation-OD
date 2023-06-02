

function [X, P] = UKFS_prediction(obj, sigma)
    % State prediction
    X = sum(obj.W(1,:).*sigma,2);

    % Covariance prediction
    [~, S] = qr([sqrt(obj.W(2,2:end)).*(sigma(:,2:end)-X) obj.Q].',0);

    if (obj.W(2,1) < 0)
        P = cholupdate(S, sqrt(abs(obj.W(2,1)))*(sigma(:,1)-X), '-');
    else
        P = cholupdate(S, sqrt(abs(obj.W(2,1)))*(sigma(:,1)-X), '+');
    end
end