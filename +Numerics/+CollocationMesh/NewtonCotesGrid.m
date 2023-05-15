%% Project: SBOPT %%
% Date: 05/05/23

%% Normal quadrature object %%
% Definition of the gaussian quadrature rule and collocation grid as an
% object

classdef NewtonCotesGrid < CollocationMesh.AbstractGrid
    
    methods 
        % Constructor 
        function [obj] = NewtonCotesGrid(m)
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
            N = obj.NumberPoints;
            w = zeros(N+1,1);   
            dumb = linspace(0,N,1000);
        
            % Compute the Newton-Cotes weights 
            for i = 1:N+1
                % Compute the functional 
                func = ones(1,length(dumb));
                for j = 0:N 
                    if ((i-1) ~= j)
                        func = func .* (dumb-j);
                    end
                end
        
                w(i) = (-1)^(N-(i-1))/(factorial(N-(i-1))*factorial((i-1)))*trapz(dumb,func);
            end

            obj.W = w.';
        end

        function [t, dt] = Domain(obj, t0, tf, tau)
            t  = (tf-t0) * obj.J * tau; 
            dt = (tf-t0) * obj.J * ones(1,length(tau));
        end
    end   
end