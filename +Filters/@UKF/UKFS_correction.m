
function [X, P, Sy] = UKFS_correction(obj, sigma, X, P, y, Y, z)
    % Covariance computation (uscented transform)
    aux_meas = y-Y;
    res = sqrt( abs(obj.W(2,:)) ) .* aux_meas;
    [~, Sy] = qr( [res(:,2:end) obj.R].', 0);

    if (obj.W(2,1) < 0)
        Sy = cholupdate(Sy, res(:,1), '-');
    else
        Sy = cholupdate(Sy, res(:,1), '+');
    end

    % Cross covariance
    Pxy = ( obj.W(2,:) .* (sigma-X) ) * aux_meas.';

    % Kalman update
    K = Pxy/Sy/Sy.';

    % Update
    X = X+K*(z-Y);
    U = K*Sy.';
    
    for i = 1:size(U,2)
        P = cholupdate(P, U(:,i), '-');  
    end
end