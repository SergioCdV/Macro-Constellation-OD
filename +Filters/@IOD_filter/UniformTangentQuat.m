
function [samples] = UniformTangentQuat(obj, L, M, mode)
    % Preallocation
    samples = zeros(4, L * M + 1);          % Output

    % Sanity checks 
    if (~exist('L', 'var'))
        L = obj.L;
    end

    if (~exist('M', 'var'))
        M = obj.M;
        beta = [obj.grid; zeros(1,M)];
    else
        beta = [obj.UniformSphere(M); zeros(1,M)];    
    end

    % Distribution on the r-sphere S^3
    One = [0,0,0,1].'; 

    % Main computation
    for i = 1:L
        samples(:,1 + M * (i-1):M*i) = (i*pi/(2*L)) * beta;
    end

    for i = 1:size(samples,2)-1
        samples(:,i) = QuaternionAlgebra.left_isoclinic(mode) * QuaternionAlgebra.exp_map(samples(:,i), One);
    end

    samples(:,end) = mode;
end