
function [State, Sigma, Pyy] = UKFA_correction(obj, sigma, State, Sigma, Y, y, z)
    % Update
    if (~any(isnan(y)))
        % Covariances matrices        
        Pxy = 0; 
        Pyy = obj.R; 
        for i = 1:size(sigma,2)
            dumb = Y(:,i)-y;
            Pxy = Pxy + obj.W(2,i) * (sigma(:,i)-State) * dumb.';
            Pyy = Pyy + obj.W(2,i) * (dumb * dumb.');
        end
        
        % Kalman gain
        K = Pxy*(Pyy\eye(size(Pyy)));

        State = State+K*(z-y);     
        Sigma = Sigma-Pxy*K.'-K*Pxy.'+K*Pyy*K.';    % Joseph update 
    else
        Pyy = [];
    end
end