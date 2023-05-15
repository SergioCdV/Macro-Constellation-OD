%% Project: SBOPT %%
% Date: 05/05/23

%% Chebyshev quadrature object %%
% Definition of the Chebyshev quadrature rule and collocation grid as an
% object

classdef ChebyshevGrid < Numerics.CollocationMesh.AbstractGrid
    
    methods 
        % Constructor 
        function [obj] = ChebyshevGrid(m)
            if (exist('m', 'var'))
                obj.NumberPoints = m;
            else
                obj.NumberPoints = 100;
            end

            if (mod(obj.NumberPoints,2) ~= 0)
                obj.NumberPoints = obj.NumberPoints + 1;
            end
 
            % Generation of the grid and grid rules
            [obj] = obj.CollocationGrid();           % Collocation grid
            [obj] = obj.QuadWeights();               % Quadrature weigths
            [obj] = obj.DiffMatrix();                % Differentiation matrix
            obj.J = 0.5;                             % Jacobian domain transformation 
        end

        % Particular methods
        function [obj] = CollocationGrid(obj)
            % Chebyshev nodes 
            i = obj.NumberPoints:-1:0;
            obj.tau = cos(pi*i/obj.NumberPoints);
        end

        function [obj] = DiffMatrix(obj)
            % Order of the matrix
            N = length(obj.tau)-1;
            
            % Main computation
            if (N == 0) 
                obj.D = zeros(N);
            else
                x = cos(pi*(0:N)/N)';
                c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
                X = repmat(x,1,N+1);
                dX = X-X';
                obj.D = (c*(1./c)')./(dX+(eye(N+1)));       % Off-diagonal entries
                obj.D = obj.D - diag(sum(obj.D,2));         % Diagonal entries
                obj.D = -obj.D;                             % Domain flip
            end  
        end

        function [obj] = QuadWeights(obj)
           % Constants 
           N = obj.NumberPoints+1;
           n = N-1; 

           % FFT preallocation of the CC weights and nodes
           c = zeros(N,2);
           c(1:2:N,1) = (2./[1 1-(2:2:n).^2]).';
           c(2,2) = 1; 
           f = real(ifft([c(1:N,:); c(n:-1:2,:)]));
           w = [f(1,1); 2*f(2:n,1); f(N,1)];
           obj.tau = n*f(1:N,2).';

           % Quadrature
           [obj.tau,i] = sort(obj.tau); 
           obj.W = w(i).';
        end

        function [t, dt] = Domain(obj, t0, tf, tau)
            t  = (tf - t0) * obj.J * (1+tau); 
            dt = (tf - t0) * obj.J * ones(1,length(tau));
        end
    end   
end