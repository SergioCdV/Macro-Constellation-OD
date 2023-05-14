
function [q] = SteepestQuat(obj, q0, w, dq)
    % Set up the Newton Rhapson method 
    maxIter = 100; 
    iter = 1; 
    tol = 1e-3; 
    GoOn = true; 
    q = q0;

    while (GoOn && iter < maxIter)
        % Compute the update term 
        dQ = zeros(4,1);
        for i = 1:size(dq,2)
            dQ = dQ + w(i) * QuaternionAlgebra.log_map(dq(:,i),q);
        end

        qn = QuaternionAlgebra.exp_map(dQ,q);
        res = q.' * qn;
        q = qn;

        % Convergence analysis
        if (abs(1-res) < tol)
            GoOn = false;
        else
            iter = iter + 1;
        end
    end
end