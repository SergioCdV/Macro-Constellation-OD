%% Project: SBOPT %%
% Date: 05/05/23

%% Bezier quadrature object %%
% Definition of the Bezier quadrature rule and collocation grid as an
% object

classdef OBezierGrid < CollocationMesh.AbstractGrid
    
    methods 
        % Constructor 
        function [obj] = OBezierGrid(m)
            if (exist('m', 'var'))
                obj.NumberPoints = m;
            else
                obj.NumberPoints = 100;
            end
 
            % Generation of the grid and grid rules
            [obj] = obj.CollocationGrid();           % Collocation grid
            [obj] = obj.QuadWeights();               % Quadrature weigths
            [obj] = obj.DiffMatrix();                % Differentiation matrix

            % Re-scaling
            obj.W = obj.W/2; 
            obj.tau = 0.5*(obj.tau+1);
            obj.J = 1;                               % Jacobian domain transformation 
        end

        % Particular methods
        function [obj] = CollocationGrid(obj)
            N = obj.NumberPoints;
            beta = .5./sqrt(1-(2*(1:N)).^(-2));         % 3-term recurrence coeffs
            T = diag(beta,1) + diag(beta,-1);           % Jacobi matrix
            [~,D] = eig(T);                             % Eigenvalue decomposition
            x = diag(D);                                % Eigenvalue matrix  
            obj.tau = sort(x);                          % Legendre nodes
        end

        function [obj] = DiffMatrix(obj)
            % Constants 
            N = obj.NumberPoints;
            L = PolynomialBases.Legendre().basis(N, obj.tau);
                        
            % Differentiation matrix
            D = zeros(N+1);
            for i = 0:N 
                for j = 0:N
                    if (i ~= j)
                        D(i+1,j+1) = (L(end,i+1)/L(end,j+1))/(obj.tau(i+1)-obj.tau(j+1));
                    elseif (i == j && j == N)
                        D(i+1,j+1) = (N*(N+1)/4);
                    elseif (i == j && j == 0)
                        D(i+1,j+1) = -(N*(N+1)/4);
                    else 
                        D(i+1,j+1) = 0; 
                    end
                end
            end

            obj.D = D;
        end

        function [obj] = QuadWeights(obj)
            % Truncation + 1
            N = obj.NumberPoints;
            N1 = N + 1;
        
            % CGL nodes
            tau = cos(pi*(0:N)/N)';
        
            % Preallocation of the Legendre Vandermonde Matrix
            P = zeros(N1);
        
            % Compute P_(N) using the recursion relation. Compute its first and second derivatives and  update x using the Newton-Raphson method.
            xold = 2;
        
            while (max(abs(tau-xold)) > eps)
                % Initialization
                xold = tau;
                P(:,1) = 1;    
                P(:,2) = tau;
                
                for k = 2:N
                    P(:,k+1) = ( (2*k-1)*tau.*P(:,k)-(k-1)*P(:,k-1) )/k;
                end
                 
                tau = xold-( tau.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
            end
        
            % Legendre basis 
            obj.tau = flip(tau).';
            L = PolynomialBases.Legendre().basis(N, tau.');
        
            % Quadrature weights
            obj.W = 2./(N*N1*L(end,:).^2);
        end

        function [t, dt] = Domain(obj, t0, tf, tau)
            t  = (tf - t0) * obj.J * tau; 
            dt = (tf - t0) * obj.J * ones(1,length(tau));
        end
    end  
end