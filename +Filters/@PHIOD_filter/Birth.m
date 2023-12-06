function [born_particles] = Birth(obj)
    % Constants 
    J = obj.M;          % Number of birth particles
    pos = 8;            % State dimensions

    % Preallocation 
    born_particles = zeros(pos + (pos-1)^2 + 1, J);

    % Uniform orbital planes birth
    born_particles(1:4,:) = QuaternionAlgebra.UniformQuat(J);

    % Action set 
    box(1,:) = [obj.Lmin obj.Lmax];
    box(2,:) = [sqrt(1-obj.emax^2) 1];
    box(3,:) = [-1 1];

    L = box(1,1) * ones(1,J) + diff(box(1,:)) .* rand(1,J);
    eta = box(2,1) * ones(1,J) + diff(box(2,:)) .* rand(1,J);
    
    born_particles(5,:) = L; 
    born_particles(6,:) = L .* eta;

    for i = 1:size(born_particles,2)
        R = QuaternionAlgebra.Quat2Matrix(born_particles(1:4,i));
        born_particles(7,i) = born_particles(6,i) .* R(end,end); 
    end

    % Anomaly distribution initialization
    born_particles(8,:) = linspace(0, 2*pi, J);

    % Covariance matrix
    sigma = blkdiag(1e-3 * eye(pos-2), deg2rad(2));
    sigma = reshape(sigma, [], 1);

    % Particles as wrapped Gaussian kernels 
    born_particles(pos+1:end-1,:) = repmat(sigma, 1, J);

    % Associated weights 
    born_particles(end,:) = repmat(1/J, 1, J);
end