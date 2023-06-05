% Compute the wrapped normal distribution
function [f] = wrapped_normal(obj, error_tol, M, mu, sigma)
    % Compute the error bound 
    n(1) = max(1+sqrt(-log(4*pi^3*error_tol^2)*sigma),1+sqrt(sigma/2)/pi);
    n(2) = max(sqrt(-log(2*pi^2*sigma*error_tol^2)/sigma),sqrt(2)/pi);
    n = ceil(n);
    
    f = zeros(length(M),1);

    % Compute the wrapped normal distributions
    if (min(n) == n(1))
        % Compute the wrapped normal
        for i = -n:n
            f = f+exp(-0.5*(M-mu+2*pi*i).^2/sigma);
        end

        f = f/sqrt(sigma*2*pi);
    else
        rho = exp(-sigma/2);
        for i = 1:n
            f = f+rho^(i^2)*cos(i*(M-mu));
        end
        f = (1+2*f)/(2*pi);
    end
end