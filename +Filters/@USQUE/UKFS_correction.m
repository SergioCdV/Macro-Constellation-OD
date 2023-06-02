
function [State, Sigma, Sy] = UKFS_correction(obj, sigma, State, Sigma, y, Y, z)
    % Covariance computation
    Sy = [sqrt(obj.W(2,2:end)).*(Y(:,2:end)-y) obj.R];
    [~, Sy] = qr(Sy.', 0);

    if (obj.W(2,1) < 0)
        Sy = cholupdate(Sy, sqrt(abs(obj.W(2,1)))*(Y(:,1)-y), '-');
    else
        Sy = cholupdate(Sy, sqrt(abs(obj.W(2,1)))*(Y(:,1)-y), '+');
    end

    Pxy = (sigma-State)*diag(obj.W(2,:))*(Y-y).';

    % Kalman update
    K = Pxy/Sy/Sy.';

    % Update
    State = State+K*(z-y);
    U = K*Sy.';
    Sigma = cholupdate(Sigma, U(:,1), "-");  
end