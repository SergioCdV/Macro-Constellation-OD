

function [particles, weights] = Pruning(obj, particles, weights)
    Set = weights > obj.PruneThresh;
    MRP = [QuaternionAlgebra.MPR2Quat(obj.a, obj.f, particles(1:4,:), false); particles(5:end,:)];

    if (any(Set))
        l = 0;
        while (any(Set))
            l = l+1; 
            [~, index] = sort(weights(Set));
            
            max = MRP(1:7,Set);
            res = max-max(:,index(end));
            P = reshape( MRP(8:end,Set), [7 7 * size(max,2)] );

            Mergeable = zeros(1, size(max,2));
            for i = 1:size(max,2)
                Mergeable(i) = res(:,i).' * (P(:, 1 + 7*(i-1): 7*i) \ res(:,i)) <= obj.MergeThresh;
            end
            Mergeable = logical(Mergeable);

            aux = weights(Set);
            weights(l) = sum( aux(Mergeable) );
            
            meanMRP = sum( aux(Mergeable) / weights(l) .* max(:, Mergeable), 2);
            
            Pm = zeros(7,7);
            for i = 1:size(max,2)
                if (Mergeable(i))
                    cov_index = 1 + 7 * (i-1) : 7 * i;
                    diff = meanMRP-max(:,i);
                    Pm = Pm + aux(i) / weights(l) * ( P(:, cov_index) + diff * diff.' );
                end
            end

            MRP(1:7,l) = meanMRP;
            MRP(8:end,l) = reshape(Pm, [], 1);
            
            q = 1;
            for k = 1:length(Set)
                if (Set(k))
                    Set(k) = ~Mergeable(q);
                    q = q+1;
                end
            end
        end
    end

    particles = [QuaternionAlgebra.MPR2Quat(obj.a, obj.f, MRP(1:3,1:l), true); MRP(4:end,1:l)];
    weights = weights(1:l);

    % Check if there are more than Jmax components 
    if (length(weights) > obj.Jmax)
        [~,index] = sort(weights); 
        index = index(end-obj.Jmax+1:end);
        weights = weights(:,index);
        particles = particles(:,index);
    end            
end