

classdef WahbaSolver
    properties
    end

    methods 
        function [obj] = WahbaSolver()
        end

        % Solvers
        [q, Sigma] = Davenports(weights, b, r);     % Davenport's q method 
        [q, Sigma] = QUEST(weights, b, r);          % Quest method
        [q, Sigma] = ESOQ(weights, b, r);           % Estimation of the optimal quaternion
        [q, Sigma] = ESOQ2(weights, b, r);          % Second estimation of the optimal quaternion
        [q, Sigma] = TRIAD(weights, b, r);          % TRIAD

    end
end