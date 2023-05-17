

function [particles, weights] = Pruning(obj, particles, weights)
    Set = weights > obj.PruneThresh;
    if (any(Set))
        l = 0;
        while (any(Set))
            l = l+1; 
            [~, index] = sort(weights(Set));
            max = mod(particles(1,Set), 2*pi);
            P = particles(2,Set);
            Mergeable = (max-max(index(end))).^2./P <= obj.MergeThresh;
            aux = weights(Set);
            w = sum(aux(Mergeable));

            particles(1,l) = dot(aux(Mergeable),max(Mergeable))/w(l);
            particles(2,l) = dot(aux(Mergeable),(P(Mergeable)+(particles(1,l)-max(Mergeable)).^2))/w(l);
            
            Set = Set(~Mergeable);
        end
    end

    % Check if there are more than Jmax components 
    if (length(weights) > obj.Jmax)
        [~,index] = sort(weights); 
        index = index(end-obj.Jmax+1:end);
        weights = weights(:,index);
        particles(1,:) = particles(1,index);
        particles(2,:) = particles(2,index);
    end
end