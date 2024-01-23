

function [State, Sigma] = PropagationStep(obj, time_step)
    % Propagation of sigma points 
    [State, Sigma] = obj.StateModel(obj.Q, obj.State, obj.Sigma, time_step);
end