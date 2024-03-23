

% Function to handle orbital elements singularities when converting from the inertial state vector
function [elements] = rv_singularity(ev, n, r, Q, elements)
    % State variables 
    e = elements(2,:);           % Orbit eccentricity
    i = elements(4,:);           % Orbit inclination 
    M = elements(6,:);           % Mean anomaly
    tol = 1E-10;                 % Circular orbit tolerance

    % Module of the position vector
    R = sqrt( dot(r(1:3,:), r(1:3,:), 1) );    
    
    % Singularity warnings 
    for j = 1:size(r,2)
        if ( any( isnan( Q(:,1+3*(j-1):3*j) ) ) )
            warning('Euler angles are numerically ill-conditioned');
        end
    end
    
    if (any(abs(e) < tol))
        warning('Orbit is circular to numerical precision');
    end
    
    if (any(abs(i) < tol))
        warning('Orbit is equatorial to numerical precision');
    end
    
    % Circular equatorial orbit singularity
    idx = abs(e) < tol && abs(i) < tol && r(2,:) >= 0; 
    M(idx) = acos( r(1,idx) ./ R(idx) ); 

    idx = abs(e) < tol && abs(i) < tol && r(2,:) < 0; 
    M(idx) = 2*pi - acos( r(1,idx) ./ R(idx) );

    % Circular inclined orbit singularity
    idx = abs(e) < tol && abs(i) >= tol && r(3,:) >= 0;
    M(idx) = acos( dot(n, r, 1) / (sqrt(dot(n, n, 1)) .* R) ); 

    idx = abs(e) < tol && abs(i) >= tol && r(3,:) < 0;
    M(idx) = 2*pi - acos( dot(n(:,idx), r(:,idx), 1) / (sqrt(dot(n(:,idx), n(:,idx), 1)) .* R(:,idx)) );

    % Equatorial orbit 
    idx = abs(e) >= tol && abs(i) < tol && ev(2,:) >= 0;
    M(idx) = acos( ev(1,idx) ./ sqrt(dot(ev(:,idx), ev(:,idx), 1)) );

    idx = abs(e) >= tol && abs(i) < tol && ev(2,:) < 0;
    M(idx) = 2*pi - acos( ev(1,idx) ./ sqrt(dot(ev(:,idx), ev(:,idx), 1)) );
    
    % Reconstruction of the orbital elements 
    elements(6,:) = M;
end