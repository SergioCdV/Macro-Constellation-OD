

function [s] = Delaunay2COE(mu, x, direction)

    if (direction)
        % Delaunay set to the COE state transformation
        M = x(1);           % Mean anomaly
        omega = x(2);       % AoP
        Omega = x(3);       % RAAN
        L = x(4);           % Delaunay action
        G = x(5);           % Angular momentum
        H = x(end);         % Polar angular momentum
    
        % Keplerian functions
        e = real(sqrt(1-(G/L)^2));
        i = acos(H/G); 
        a = L^2/mu;
        
        s = real([a e Omega i omega M]);

    else
        % COE set to the Delaunay state transformation
        a = x(1);       % Semimajor axis
        e = x(2);       % Eccentricity
        Omega = x(3);   % RAAN
        i = x(4);       % Inclination
        omega = x(5);   % AoP
        M = x(6);       % Mean anomaly
    
        % Keplerian functions
        L = sqrt(mu * a);
        G = L * sqrt(1-e^2);
        H = G * cos(i); 
        
        s = real([M omega Omega L G H]);
    end
end