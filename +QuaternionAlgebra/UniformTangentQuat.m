
function [quat] = UniformTangentQuat(L, M, mode)
    % Preallocation
    quat = zeros(size(mode,1), (L * M + 1) * 1);          % Output
    One = [0;0;0;1];

    beta = [QuaternionAlgebra.UniformSphere(M); zeros(1,M)];    

    % Main computation
    for i = 1:L
        tangent(:,1 + M * (i-1):M*i) = (i*pi/(2*L)) * beta;
    end
    
    % Copy the mode
    J = L*M;
    quat(:,1) = mode;

    % Perturb the attitude
    for i = 1:J
        % QuaternionAlgebra.exp_map(tangent(:,j), mode(1:4,i));
        quat(1:4,1+i) = QuaternionAlgebra.right_isoclinic(mode) * QuaternionAlgebra.exp_map(tangent(:,i), One);
    end

    quat(1:4,:) = quat(1:4,:) ./ sqrt(dot(quat(1:4,:), quat(1:4,:),1));
end