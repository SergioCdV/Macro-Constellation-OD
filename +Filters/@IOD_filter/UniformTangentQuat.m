
function [samples] = UniformTangentQuat(obj, L, M, mode)
    % Preallocation
    samples = zeros(4, L * M + 1);          % Output

    % Sanity checks 
    if (~exist('L', 'var'))
        L = obj.L;
    end

    % Distribution on the r-sphere S^3
    if (~exist('M', 'var'))
        M = obj.M;
        beta = [obj.grid; zeros(1,M)];
    else
        beta = [obj.UniformSphere(M); zeros(1,M)];    
    end

    % Main computation
    for i = 1:L
        tangent(:,1 + M * (i-1):M*i) = (i*pi/(2*L)) * beta;
    end

    for i = 1:size(samples,2)-1
        % samples(:,i) =  QuaternionAlgebra.right_isoclinic(mode) * QuaternionAlgebra.exp_map(tangent(:,i), ones);
        samples(:,i) =  QuaternionAlgebra.exp_map(tangent(:,i), mode);
    end

    samples(:,end) = mode;
    samples = samples ./ sqrt(dot(samples(1:4,:), samples(1:4,:),1));
end