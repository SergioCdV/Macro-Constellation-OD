

function [particles, weights] = Initialization(obj)
    % Sanity checks 
    if (isempty(obj.N) || isempty(obj.L) || isempty(obj.M))
        error('The internal filter configuration has not been completed.');
    else  
        % Generate the initial multi-modal uniform distribution on the quaternion sphere 
        quat = obj.UniformQuat(obj.N);
    
        % Generate the initial uniform distribution on the action space
        mu = [1.08; 0.001; cos(deg2rad(50))];
        sigma = diag([1 0.001 1]);
        sigma = zeros(3);
    
        % Generate the particles
        particles = zeros(7, obj.N * (obj.L * obj.M + 1));
        for i = 1:obj.N
            % Generate the actions
            actions = mvnrnd(mu, sigma, (obj.L * obj.M + 1)).';

            % Physical constraint on the Delaunay action
            actions(1,:) = repmat(mu(1), 1, size(actions,2)) .* (actions(1,:) <= 1) + actions(1,:) .* (actions(1,:) > 1);
    
            % Physical constraint on the eccentricity
            actions(2,:) = repmat(mu(2), 1, size(actions,2)) .* (actions(2,:) < 0) + actions(2,:) .* (actions(1,:) >= 0);
    
            % Correlation
            actions(2,:) = actions(1,:) .* sqrt(1-actions(2,:).^2);
            actions(3,:) = actions(2,:) .* actions(3,:);

            % particles(:,1+(obj.L * obj.M + 1)*(i-1):(obj.L * obj.M + 1)*i) = [obj.UniformTangentQuat(obj.L, obj.M, quat(1:4,i)); actions];
            particles(5:7,1+(obj.L * obj.M + 1)*(i-1):(obj.L * obj.M + 1)*i) = actions;
        end
        particles(1:4,:) = obj.UniformQuat(obj.N * (obj.L * obj.M + 1));

        % Uniform generation of the weights 
        weights = repmat(1/size(particles,2), 1, size(particles,2));
    end
end