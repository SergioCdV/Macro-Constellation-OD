

function [particles, weights] = Pruning(obj, particles, weights)
    Set = weights > obj.PruneThresh;
    if (any(Set))
        l = 0;
        while (any(Set))
            l = l+1; 
            [~, index] = sort(weights(Set));
            max = particles(2:8,Set);
            P = particles(9:end,Set);
            res = max-max(:,index(end));
            Mergeable = res.' * res ./ reshape(P, [7 7]) <= obj.MergeThresh;
            aux = weights(Set);
            w(l) = sum(aux(Mergeable));

            particles(:,l) = particles(:,index(end));
            particles(8,l) = dot(aux(Mergeable),max(Mergeable))/w(l);
            particles(end,l) = dot(aux(Mergeable),(P(Mergeable)+(particles(1,l)-max(Mergeable)).^2))/w(l);
            
            q = 1;
            for k = 1:length(Set)
                if (Set(k))
                    Set(k) = ~Mergeable(q);
                    q = q+1;
                end
            end
        end
    end

    % Check if there are more than Jmax components 
    if (length(weights) > obj.Jmax)
        [~,index] = sort(weights); 
        index = index(end-obj.Jmax+1:end);
        weights = weights(:,index);
        particles = particles(:,index);
    end
end