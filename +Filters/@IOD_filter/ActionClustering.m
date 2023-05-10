

function [c, Sigma] = ActionClustering(obj, samples, N)
    % Configuration 
    rng(1); 

    % Perform K-means clustering over the action data
    [index, c] = kmeans(samples.', N);

    % Compute the standard deviation 
    Sigma = zeros(N,1);
    for i = 1:N
        sigma = (samples(:,index(index == i))-c(i)).^2;
        Sigma(i,1) = sum(sigma)/sum(index(index == i));
    end
end