%% Project: SBOPT %%
% Date: 05/05/23

%% Trapezoidal quadrature object %%
% Definition of the gaussian quadrature rule and collocation grid as an
% object

classdef TrapezoidalGrid < CollocationMesh.AbstractGrid
    
    methods 
        % Constructor 
        function [obj] = TrapezoidalGrid(m)
            if (exist('m', 'var'))
                obj.NumberPoints = m;
            else
                obj.NumberPoints = 100;
            end
 
            % Generation of the grid and grid rules
            [obj] = obj.CollocationGrid();           % Collocation grid
            [obj] = obj.QuadWeights();               % Quadrature weigths
            [obj] = obj.DiffMatrix();                % Differentiation matrix
            obj.J = 1;                               % Jacobian domain transformation 
        end

        % Particular methods
        function [obj] = CollocationGrid(obj)
            obj.tau = linspace(0, 1, obj.NumberPoints + 1);
        end

        function [obj] = DiffMatrix(obj)
            obj.D = [];
        end

        function [obj] = QuadWeights(obj)
            obj.W = [];
        end

        function [t, dt] = Domain(obj, t0, tf, tau)
            t  = (tf-t0) * obj.J * tau; 
            dt = (tf-t0) * obj.J * ones(1,length(tau));
        end
    end   
end