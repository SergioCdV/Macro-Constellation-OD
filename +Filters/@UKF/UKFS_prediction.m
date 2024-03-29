

function [X, S] = UKFS_prediction(obj, sigma)
    % State prediction
    X = sum(obj.W(1,:).*sigma,2);

    % Covariance prediction
    res = (sigma-X) * diag( sqrt(abs(obj.W(2,:))) );
    [~, S] = qr([res(:,2:end) obj.Q].',0);

    if (obj.W(2,1) < 0)
        S = cholupdate(S, res(:,1), '-');
    else
        S = cholupdate(S, res(:,1), '+');
    end
end