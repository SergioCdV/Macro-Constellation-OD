
function [samples] = UniformQuat(m)
    % Initialization 
    u = rand(3,m);

    % Random sampling of the 3 sphere 
    samples = [sqrt(u(1,:)) .* cos(2*pi*u(3,:));   ...
               sqrt(u(1,:)) .* sin(2*pi*u(3,:));   ...
               sqrt(1-u(1,:)) .* sin(2*pi*u(2,:)); ...
               sqrt(1-u(1,:)) .* cos(2*pi*u(2,:))];
end