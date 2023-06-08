function [X] = StateEstimation(obj, samples, weights, T)
    % Configuration 
    rng(1); 

    % Extract the relevant particles 
    pruned_samples = samples;
    pruned_weights = weights;

    % Perform K-means clustering over the quaternions and perform the state estimation for the perifocal attitude
    [c, ~, index] = obj.QuatClustering(pruned_weights, pruned_samples(1:4,:));

    % Estimate the number of planes based on the quaternion distribution
    N = size(c,2);

    % Preallocation 
    X = zeros(7+6^2, N);
    pos = 0;

    for i = 1:max(index)
        ID = index == i;
        pos = i;

        % Quaternion estimation
        [X(1:4,i), Sigma_q] = QuaternionAlgebra.AverageQuat(ones(1,sum(ID)), pruned_samples(1:4,ID));

        % Estimation of the action set
        w = pruned_weights(1,ID) / sum(pruned_weights(1,ID));
        X(5:7,i) = sum(w .* pruned_samples(5:7,ID), 2);

        % Covariance
        qs = pruned_samples(1:4,ID);
        a = zeros(3,size(qs,2));
        Q = QuaternionAlgebra.right_isoclinic(X(1:4,i));
        for j = 1:size(qs,2)
            dq = Q * QuaternionAlgebra.quaternion_inverse(qs(1:4,j));
            a(:,j) = dq(1:3,1) / (1 + dq(4,1));
        end

        mu = [zeros(3,1); X(5:7,i)];
        s = [a; pruned_samples(5:7,ID)];
 
        sigma = (w.*s-mu*w) * (w.*s-mu*w).';

        if (sum(ID) > 1)
            sigma = sum(pruned_weights(1,ID),2) / (sum(pruned_weights(1,ID),2)^2-sum(pruned_weights(1,ID).^2,2)) * sigma;
        else
            sigma = zeros * eye(6);
        end

        sigma(1:3,1:3) = Sigma_q;
        sigma = 0.5 * (sigma + sigma.') + obj.PD_tol * eye(6);
        X(8:end,i) = reshape(sigma, [], 1);
    end

    if (any(index < 0))
        ID = index == -1;
        
        for i = pos+1:pos+sum(ID)
            % Quaternion estimation
            [X(1:4,i), Sigma_q] = QuaternionAlgebra.AverageQuat(ones(1,sum(ID)), pruned_samples(1:4,ID));
    
            % Estimation of the action set
            w = pruned_weights(1,ID) / sum(pruned_weights(1,ID));
            X(5:7,i) = sum(w .* pruned_samples(5:7,ID), 2);
    
            % Covariance 
            Q = QuaternionAlgebra.quaternion_inverse(X(1:4,i));
            qs = pruned_samples(1:4,ID);
            a = zeros(3,size(qs,2));
            for j = 1:size(qs,2)
                dq = QuaternionAlgebra.right_isoclinic(qs(1:4,j)) * Q;
                aux = QuaternionAlgebra.log_map(dq, [0;0;0;1]);
                a(:,j) = aux(1:3);
            end
    
            mu = [zeros(3,1); X(5:7,i)];
            s = [a; pruned_samples(5:7,ID)];
    
            sigma = (w.*s-mu*w) * (w.*s-mu*w).';
    
            if (sum(ID) > 1)
                sigma = sum(pruned_weights(1,ID),2) / (sum(pruned_weights(1,ID),2)^2-sum(pruned_weights(1,ID).^2,2)) * sigma;
            else
                sigma = zeros * eye(6);
            end
    
            sigma(1:3,1:3) = Sigma_q;
            sigma = 0.5 * (sigma + sigma.') + obj.PD_tol * eye(6);
            X(8:end,i) = reshape(sigma, [], 1);
        end
    end
end
