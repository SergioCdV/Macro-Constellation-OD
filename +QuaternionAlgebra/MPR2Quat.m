
function [X] = MPR2Quat(a, f, x, direction)

    if (direction)
        aux = dot(x(1:3,:), x(1:3,:), 1);
        q = (-a * aux + f * sqrt(f^2 + (1-a^2) * aux) ) ./ (f^2 + aux );
        q = [(a + q(1,:)) / f .* x(1:3,:); q(1,:)];
        X = q;
    else
        X = f * x(1:3,:) ./ (a + x(4,:));
    end
end