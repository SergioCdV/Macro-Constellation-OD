

function [State, Sigma] = CorrectionStep(obj, State, Sigma)        
    % Measurement prediction 
    [y, H] = obj.ObservationModel(State);

    % Kalman gain 
    P = (obj.R+H*Sigma*H.');
    K = Sigma*H.'*P^(-1);

    % Update
    State = State+K*(z-y);
    Sigma = (eye(obj.StateDim)-K*H)*Sigma*(eye(obj.StateDim)-K*H).'+K*obj.R*K.';
end