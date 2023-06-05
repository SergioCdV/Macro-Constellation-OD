
function [State, Sigma, Sy] = UKFS_correction(obj, sigma, State, Sigma, Y, y, z)
    % Covariance computation
    res = (Y-y) * diag( sqrt(abs(obj.W(2,:))) );
    [~, Sy] = qr([res(:,2:end) obj.R].',0);

    if (obj.W(2,1) < 0)
        Sy = cholupdate(Sy, res(:,1), '-');
    else
        Sy = cholupdate(Sy, res(:,1), '+');
    end

    Pxy = 0; 
    for i = 1:size(sigma,2)
        Pxy = Pxy + obj.W(2,i) * (sigma(:,i)-State) * (Y(:,i)-y).';
    end

    % Kalman update
    K = Pxy/Sy/Sy.';

    % Update
    State = State+K*(z-y);
    U = K*Sy.';
    
    for i = 1:size(U,2)
        Sigma = cholupdate(Sigma, U(:,i), "-");  
    end
end