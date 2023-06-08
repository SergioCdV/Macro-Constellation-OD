
function [samples] = AffineSampling(obj, m, mu, Sigma)

    % Cholesky factorization of Sigma 
    [A, flag] = chol(Sigma);
    M = size(mu,1);
    if (~flag)         
        if (mod(M,2) ~= 0)
            M = M+1;

            % Generate the uniform numbers 
            u = rand(M, m);

            % Transform them to univariate normal distributions
            for i = 1:size(mu,1)-1
                u(i:i+1,:) = sqrt(-2 * log(u(i,:))) .* [cos(2*pi*u(i+1,:)); sin(2*pi*u(i+1,:))];
            end

            u = u(1:end-1,:);

        else
            % Generate the uniform numbers 
            u = rand(M, m);

            % Transform them to univariate normal distributions
            for i = 1:size(mu,1)-1
                u(i:i+1,:) = sqrt(-2 * log(u(i,:))) .* [cos(2*pi*u(i+1,:)); sin(2*pi*u(i+1,:))];
            end
        end

        % Compute the final samples
        samples = real(mu + A.' * u);           % Attention with the regularized covariance
    else
        warning('Covariance matrix is not positive definite.');
    end
end