
function [s] = Delaunay2Lara(x, direction)
    if (direction)
        M = x(1);           % Mean anomaly
        omega = x(2);       % AoP
        Omega = x(3);       % RAAN
        L = x(4);           % Delaunay action
        G = x(5);           % Angular momentum
        H = x(end);         % Polar angular momentum

        e = real(sqrt(1-(G/L)^2));
        cos_i = H/G;

        % Solve for the true anomaly
        [nu, E] = Astrodynamics.KeplerSolver(e,M);
    
        % Compute the radial and velocity distance in the perifocal frame
        R = G^2./(1+e*cos(nu));

        theta = nu + omega;
        
        s(1,1) = Omega + theta;
        s(2,1) = sin(theta) * sqrt(1-cos_i^2);
        s(3,1) = cos(theta) * sqrt(1-cos_i^2);
        s(4,1) = R;
        s(5,1) = L/R * e * sin(E);
        s(6,1) = G;
    else 
        s = [];
    end
end