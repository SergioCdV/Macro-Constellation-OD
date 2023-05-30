
function [samples, a] = PerifocalQuatSampling(obj, particles)
    % Preallocation
    pos = 8;
    J = size(particles,2);
    samples = particles;          % Output

    % Gibbs covariance 
    Sigma = reshape(particles(pos+1:pos+(pos-1)^2,1), [pos-1 pos-1]);

    % Generate the Gibbs vector 
    a = mvnrnd(zeros(3,1), Sigma(1:3,1:3), J).';

    for i = 1:size(samples,2)
        samples(1:4,i) = QuaternionAlgebra.exp_map([a(:,i); 0], samples(1:4,i));
    end

    % Two expensive
%     J = obj.L*obj.R+1;
%     samples = zeros(size(particles,1), J * size(particles,2));    % Output
%     
%     for i = 1:size(particles,2)
%         % Copy the mode
%         samples(1:end-1, 1+J*(i-1):J*i) = repmat(particles(1:end-1,i), 1, J);
%         samples(end, 1+J*(i-1):J*i) = repmat(particles(end,i)/J, 1, J);
% 
%         % Perturb the perifocal attitude
%         samples(1:4,1+J*(i-1):J*i) = QuaternionAlgebra.UniformTangentQuat(obj.L, obj.R, particles(1:4,i));       
%     end

    % Numerical conditioning
    samples(1:4,:) = samples(1:4,:) ./ sqrt(dot(samples(1:4,:), samples(1:4,:),1));
end