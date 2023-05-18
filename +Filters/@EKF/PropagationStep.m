

function [State, Sigma] = PropagationStep(obj, time_step)
    % Propagation of sigma points 
    [State, Sigma] = obj.StateModel(obj.Q, time_step, obj.State, obj.Sigma);
end