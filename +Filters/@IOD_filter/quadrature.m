
function [Y] = quadrature(obj, x)
    % CC quadrature
    % Y = dot(obj.dt(2,:) .* obj.W, x, 2);

    % MC integration 
    Y = sum(x) * 2*pi/length(obj.nu);
end