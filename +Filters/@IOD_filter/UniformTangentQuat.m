
function [samples] = UniformTangentQuat(obj, L, M)
    % Preallocation
    samples = zeros(4, L * M);          % Output
    r = linspace(0, pi/2, L);           % Radius 

    % Main computation
    for i = 1:length(r)
        % Distribution on the r-sphere S^3
        samples(1:3, 1 + M*(i-1):M*i) = r(i) * obj.UniformSphere(M);
    end
end