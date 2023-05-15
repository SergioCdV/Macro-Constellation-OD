
classdef (Abstract) AbstractGrid
    properties
        NumberPoints;   % Number of points in the grid
        tau;            % Collocation grid
        W;              % Quadrature weights
        J;              % Domain transformation Jacobian
        D;              % Differentiation matrix
    end

    methods 
        % Class methods
        [obj] = CollocationGrid(m);
        [obj] = DiffMatrix(m);
        [obj] = QuadWeights(m);
        [t] = Domain(t0, tf, tau);
    end
end