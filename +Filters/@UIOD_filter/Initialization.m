

function [particles, weights] = Initialization(obj)
    % Sanity checks 
    if (isempty(obj.N) || isempty(obj.M))
        error('The internal filter configuration has not been completed.');
    else  
        % Number of mixture components per plane
        J = max(obj.Jmax, obj.M); 
        pos = 7;

        % Preallocation 
        particles = zeros(pos + (pos-1)^2, J);

        % Quaternion distribution initialization
        particles(1:4,:) = QuaternionAlgebra.UniformQuat(J);

        % Action set 
        box(1,:) = [obj.Lmin obj.Lmax];
        box(2,:) = [sqrt(1-obj.emax^2) 1];
        box(3,:) = [-1 1];

        L = box(1,1) * ones(1,J) + diff(box(1,:)) .* rand(1,J);
        eta = box(1,1) * ones(1,J) + diff(box(2,:)) .* rand(1,J);
        
        particles(5,:) = L; 
        particles(6,:) = L .* eta; 

        for i = 1:size(particles,2)
            R = QuaternionAlgebra.Quat2Matrix(particles(1:4,i));
            particles(7,i) = particles(6,i) .* R(end,end); 
        end

        % Covariance matrix
        sigma = 1e-3 * eye(pos-1);
        sigma = reshape(sigma, [], 1);

        % Particles as wrapped Gaussian kernels 
        particles(pos+1:end,:) = repmat(sigma, 1, J);

        % Uniform generation of the weights 
        weights = repmat(obj.N/size(particles,2), 1, size(particles,2));
    end
end