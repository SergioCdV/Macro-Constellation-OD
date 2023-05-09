
function [samples] = AffineSampling(obj, m, mu, Sigma)

    % Cholesky factorization of Sigma 
    [A, flag] = chol(Sigma);

    if (~flag)    
        % Generate the uniform numbers 
        u = rand(size(mu,1), m);

        if (mod(m,2) ~= 0)
            m = m+1;

            % Transform them to univariate normal distributions
            for i = 1:size(mu,1)-1
                u(i:i+1,:) = sqrt(-2 * log(u(i,:))) .* [cos(2*pi*u(i+1,:)); sin(2*pi*u(i+1,:))];
            end

            u = u(1:end-1,:);

        else
            % Transform them to univariate normal distributions
            for i = 1:size(mu,1)-1
                u(i:i+1,:) = sqrt(-2 * log(u(i,:))) .* [cos(2*pi*u(i+1,:)); sin(2*pi*u(i+1,:))];
            end
        end

        % Compute the final samples
        samples = mu + A.' * u;

    else 
        warning('Covariance matrix is not positive definite.');
        samples = mvnrnd(mu, Sigma, m).';
    end
end