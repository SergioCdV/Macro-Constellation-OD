

classdef WahbaSolver
    properties
    end

    methods 
        function [obj] = WahbaSolver()
        end

        % Solvers
        [q, Sigma] = Davenports(obj, weights, b, r);     % Davenport's q method 
        [q, Sigma] = QUEST(obj, weights, b, r);          % Quest method
        [q, Sigma] = ESOQ(obj, weights, b, r);           % Estimation of the optimal quaternion
        [q, Sigma] = Analytical2(obj, weights, b, r);    % Analytical solution for two measurements
        [q, Sigma] = TRIAD(obj, weights, b, r);          % TRIAD

    end
end