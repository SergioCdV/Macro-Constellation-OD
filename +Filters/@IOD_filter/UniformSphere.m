
function [samples] = UniformSphere(obj, m)
    % Initialization 
    u = mvnrnd(zeros(3,1), eye(3), m).';

    % Uniform distribution 
    samples = u./sqrt(dot(u,u,1));
end