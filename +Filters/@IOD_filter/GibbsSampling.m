
function [samples] = GibbsSampling(obj, m, mu, Sigma)  
    % Initialization 
    samples = rand(size(mu,1),m);     % Uniform sampling
    dims = 1:size(mu,1);        % Index for each dimension
 
    % Gibbs sampling
    t = 1;
    while (t < m)
        t = t + 1;
        T = [t-1,t];

        % Loop over each dimension
        for i = 1:size(mu,1) 
            nIx = dims ~= i; 

            % Conditional mean
            muCond = mu(i) + Sigma(i) * Sigma(i,nIx).' * Sigma(nIx,nIx)^(-1) * (samples(nIx,T(i))-mu(nIx));

            % Conditional covariance
            varCond = Sigma(i,i) - Sigma(i,nIx).' * Sigma(i,i)^(-1) * Sigma(i, nIx);

            % Drawing
            samples(i,t) = normrnd(muCond, varCond);
        end
    end
end