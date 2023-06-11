
function [samples, a] = PerifocalQuatSampling(obj, particles)
    % Preallocation
    pos = 8;
    J = size(particles,2);
    samples = particles;          % Output

    % Gibbs covariance 
    Sigma = reshape(particles(pos+1:pos+(pos-1)^2,1), [pos-1 pos-1]);

    % Generate the MPR vectors 
%     Sigma = Sigma(1:3,1:3) - Sigma(1:3,4:6) * Sigma(4:6,4:6)^(-1) * Sigma(1:3,4:6).';
%     [~, flag] = chol(Sigma);
%     if (flag)
%         Sigma = 0.5 * (Sigma + Sigma.') + obj.PD_tol * eye(size(Sigma)); 
%     end

    Sigma = Sigma(1:3,1:3);
    a = mvnrnd(zeros(3,1), Sigma, J).';

    for i = 1:size(samples,2)
        dq = QuaternionAlgebra.MPR2Quat(1, 1, a(:,i), true);
        samples(1:4,i) = QuaternionAlgebra.right_isoclinic( samples(1:4,i) ) * dq;
    end
    
    % Numerical conditioning
    samples(1:4,:) = samples(1:4,:) ./ sqrt(dot(samples(1:4,:), samples(1:4,:),1));
end