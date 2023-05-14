
function [samples] = GibbsSampling(obj, m, mu, Sigma, a)  
    % Initialization 
    samples = rand(size(mu,1),m);     % Uniform sampling
    dims = 1:size(mu,1);              % Index for each dimension

    if (~exist('a', 'var'))
        a = [Inf; -Inf];
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
            
            muCond = mu(i) + Sigma(i,nIx) * Sigma(nIx,nIx)^(-1) * (c-mu(nIx));

            % Conditional covariance
            varCond = Sigma(i,i) - Sigma(i, nIx) * Sigma(i,i)^(-1) * Sigma(nIx,i);

            % Drawing 
            l = (a(i,1:2) - muCond)/varCond;
            CDF = (1/2) * ( erf(l(2)) - erf(l(1)) );
            samples(i,t) = normrnd(muCond, varCond) / CDF;

            % Update the state 
            T(i) = T(i) + 1;
        end
    end
end