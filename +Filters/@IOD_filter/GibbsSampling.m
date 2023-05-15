
function [samples] = GibbsSampling(obj, m, mu, Sigma, R)  
    % Initialization 
    samples = rand(size(mu,1),m);     % Uniform sampling
    dims = 1:size(mu,1);              % Index for each dimension

    if (~exist('R', 'var'))
        R = repmat([-Inf; Inf], size(mu,1), 1);
    end
 
    % Gibbs sampling
    t = 1;
    while (t < m)
        % Update the counter 
        t = t+1;
        T = (t-1) * ones(size(mu,1),1);

        % Loop over each dimension
        for i = 1:size(mu,1) 
            nIx = dims ~= i; 

            % Conditional mean
            b = samples(nIx,:);
            c = zeros(sum(nIx),1);
            index = T(nIx);
            for j = 1:length(index)
                c(j,1) = b(j,index(j));
            end
            
            A = Sigma(nIx,nIx)^(-1);
            muCond = mu(i) + Sigma(i,nIx) * A * (c-mu(nIx));

            % Conditional covariance
            varCond = sqrt(2) * sqrt( Sigma(i,i) - Sigma(i, nIx) * A * Sigma(nIx,i) );

            % Drawing 
            u = rand();
            range = (R(i,1:2) - muCond);
            rn = range/varCond;

            if ((prod(sign(rn)) > 0) && all(abs(rn) > 4))
                varCond = sign(rn(1)) * varCond;
                rn = abs(rn);
                left = min(rn); 
                right = max(rn);
                a = 0.147;

                % Asymptotic expansion of erf and its inverse
                x2 = left*left;
                ax2 = a*x2;
                e1 = (4/pi+ax2) ./ (1+ax2);
                e1 = exp(-x2.*e1);                  % e1 < 3.0539e-008 for asymthreshold = 4
                x2 = right*right;
                ax2 = a*x2;
                e2 = (4/pi+ax2) ./ (1+ax2);
                e2 = exp(-x2.*e2);                  % e2 < 3.0539e-008 for asymthreshold = 4
                
                % Taylor series of erf(right)-erf(left) ~= sqrt(1-e2)-sqrt(1-e1)
                de = -0.5*(e2-e1) -0.125*(e2-e1)*(e2+e1);
                
                % Taylor series of erf1 := erf(left)-1 ~= sqrt(1-e1)-1 
                erf1 = (-0.5*e1 - 0.125*e1^2);
                f = erf1 + de*u;                    % = z - 1; thus z = 1+f
                l = log(-f.*(2 + f));               % log(-2f-f^2) = log(1-z.^2);
                b = 2/(pi*a) + l/2;
                X = varCond * sqrt(-b + sqrt(b.^2-l/a));
            else
                CDF = erf(range/varCond);
                X = varCond * erfinv(u * diff(CDF) + CDF(1));
            end
            
            X = max(min(X, range(2)), range(1));
            samples(i,t) = muCond + X;
            
            % Update the state 
            T(i) = T(i) + 1;
        end
    end
end