%% Constellation tracking
% Author: Sergio Cuevas
% Date: 13/03/2024

%% Kepler equation solver
% The following function solves Kepler's equation by means of the Laguerre-Conway method

% Inputs: -e, scalar, the orbital eccentricity
%         - M, scalar, the mean anomaly

% Ouput: - theta, scalar, the associated true anomaly
%        - E, scalar, the associated eccentric anomaly

function [theta, E] = KeplerSolver(e, M)
    % Laguerre-Conway's method
    maxIter = 10;       % Maximum number of iterations
    iter = 1;           % Initial iteration
    GoOn = true;        % Convergene boolean flag
    k = 5;              % Laguerre constant
    tol = 1E-15;        % Convergence tolerance

    % Warm start
    u = M + e;
    E = M .* (1-sin(u)) + u .* sin(M) ./ (1+sin(M)-sin(u));

    while (iter < maxIter && GoOn)
        fn = E - e .* sin(E) - M;
        dfn = 1 - e .* cos(E);
        ddfn = e .* sin(E);

        root = abs( (k-1)^2 * dfn^2 - k * (k-1) * fn * ddfn );
        dg(1) = dfn + sqrt(root);
        dg(2) = dfn - sqrt(root);

        dg = dg( abs(dg) == max( abs(dg) ) );

        dn = - k * fn / dg;
        E = E + dn;

        if all(abs(dn) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end

    % Solve for the true anomaly 
    sin_E = sqrt(1-e.^2) .* sin(E) ./ (1-e.*cos(E));
    cos_E = (cos(E)-e) ./ (1-e.*cos(E));
    theta = atan2(sin_E, cos_E);
end