function [theta] = KeplerSolver(e, M)
    % Laguerre-Conway's method
    maxIter = 10;
    iter = 1;
    GoOn = true;
    k = 5; 
    tol = 1e-5;

    % Warm start
    u = M + e;
    E = M*(1-sin(u)) + u*sin(M)/(1+sin(M)-sin(u));

    while (iter < maxIter && GoOn)
        fn = E -e * sin(E) - M;
        dfn = 1 - e * cos(E);
        ddfn = e * sin(E);
        dn = -k*fn/(dfn+sqrt((k-1)^2*dfn^2-k*(k-1)*dfn*ddfn));
        E = E + dn;

        if (abs(dn) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end

    % Solve for the true anomaly 
    sin_E = sqrt(1-e^2) * sin(E)/(1-e*cos(E));
    cos_E = (cos(E)-e)/(1-e*cos(E));
    theta = atan2(sin_E, cos_E);
end