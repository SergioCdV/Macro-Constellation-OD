
function [particles, weights] = Pruning(obj, particles, weights)
    if (size(particles,2) > obj.Jmax)
        Set = weights > obj.PruneThresh;
        N = sum(Set);
        if (any(Set))
            l = 0;
            while (any(Set) && l < N)
                l = l+1; 
                [~, index] = sort(weights(Set));
                max = particles(:,Set);
                Mergeable = acos(max(1:4,:).'*max(1:4,index(end))) <= obj.MergeThresh;

                if (sum(Mergeable))
                    aux = weights(Set);
                    weights(l) = sum(aux(Mergeable));
    
                    X = sum(aux(Mergeable) .*  max(1:4,Mergeable),2) / norm(sum(aux(Mergeable) .*  max(1:4,Mergeable),2));
                    particles(1:4,l) = obj.SteepestQuat(X, aux(Mergeable), max(1:4,Mergeable));
                    particles(5:7,l) = sum(aux(Mergeable) .* max(5:7,Mergeable),2) / weights(l);
                    
                    Set = Set(~Mergeable);
                end
            end
        end

        % Check if there are more than Jmax components 
        if (size(weights,2) > obj.Jmax)
            [~,index] = sort(weights); 
            index = index(end-obj.Jmax+1:end);
            weights = weights(index);
            particles = particles(:,index);
        end
    end
end