%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Cartesian state vector to classical orbital elements %%
% This file contains the function to change from the state vector to the classical orbital elements
% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - vector s, the inertial state vector to be converted (position and velocity)
%         - string frame, indicating in which reference frame the state
%           vector is expressed (inertial or perifocal)

% Ouputs: - vector elements, containing the mean classical Euler orbital elements (a, e, RAAN, i, omega, M, p)

function [elements] = ECI2COE(mu, s, direction)    
    % Main computation 
    if (direction) 
        elements = rv2coe(mu, s);
    else
        elements = coe2state(mu, s);
    end
end

%% Auxiliary functions 
% Function to compute the orbital elements from the inertial state vector 
% Inputs: - scalar mu, the gravitational parameter of the central body.
%         - vector s, containing the state vector in the inertial frame (position and velocity vectors)
% Ouputs: - vector elements, containing the classical orbital elements. 
function [elements] = rv2coe(mu, s)
    % State variables     
    r = s(1:3,:);                                   % Position vector
    v = s(4:6,:);                                   % Velocity vector 
    
    % Compute the eccentricy and angular momentum vectors 
    h = cross(r,v);                                 % Angular momentum vector
    e = cross(v,h) / mu - r ./ sqrt(dot(r,r,1));    % Eccentricity vector
    K = [0; 0; 1];                                  % Inertial Z axis unit vector
    n = cross(repmat(K, 1, size(h,2)), h);          % Node vector
    e_norm = sqrt(dot(e,e,1));                      % Norm of the eccentricity function

    % Compute orbital energy 
    H = dot(v,v,1) / 2 - mu ./ sqrt(dot(r,r,1));

    a = zeros(1,size(s,2));             % Semimajor axis
    p = zeros(1,size(s,2));             % Semilatus rectum
    
    % Determine type of orbit 
    a(1, e_norm ~= 1) = -mu ./ (2*H(1, e_norm ~= 1));                                   % Semimajor axis of the orbit
    a(1, e_norm == 1) = Inf * ones(1,sum(e_norm == 1));                                 % Semimajor axis of the orbit

    p(1, e_norm ~= 1) = a(1, e_norm ~= 1) .* (1-e_norm(e_norm ~= 1).^2);                % Semilatus rectum of the orbit
    p(1, e_norm == 1) = sqrt(dot(h(:,e_norm == 1),h(:,e_norm == 1),1)) .^2 / mu;        % Semilatus rectum of the orbit

    % Compute the unit perifocal triad  
    m = e ./ e_norm; 
    m(:, e_norm == 0) = n(:, e_norm == 0);

    k = h ./ sqrt( dot(h, h, 1) ); 
    j = cross(k,m);

    Q = reshape([m; j; k], 3, []).';  

    % Rest of elements 
    RAAN = atan2(Q(3:3:end,1),-Q(3:3:end,2));             % RAAN
    omega = atan2(Q(1:3:end,3),Q(2:3:end,3));             % Argument of perigee
    I = acos(Q(3:3:end,3));                               % Inclination

    % Position in the perifocal frame 
    r0 = r;
    for i = 1:size(r,2)
        r0(:,i) = Q(1+3*(i-1):3*i,:) * r(:,i);
    end

    % Mean anomaly
    theta = atan2(r0(2,:), r0(1,:));                      % True anomaly of the orbit
    den = 1 + e_norm .* cos(theta);
    sinE = sqrt(1-e_norm.^2) .* sin(theta) ./ den;        % Sine of the eccentric anomaly
    cosE = (e_norm+cos(theta)) ./den;                     % Cosine of the eccentric anomaly
    E = atan2(sinE, cosE);                                % Eccentric anomaly
    M = E - e_norm .* sin(E);                             % Mean anomaly
        
    % Save the classical orbital elements 
    elements = [a; e_norm; RAAN.'; I.'; omega.'; M; p];

    % Non-singular COE
    elements = Astrodynamics.rv_singularity(e, n, r, Q, elements);     
end

% Transform COE to Cartesian state vector
% Inputs: - scalar mu, the gravitational parameter of the central body.
%         - vector elements, containing the classical orbital elements. 
% Ouputs: - vector s, containing the state vector in the inertial frame (position and velocity vectors).
function [s] = coe2state(mu, elements)
    % Constants 
    e = elements(2,:);           % Eccentricity of the orbit
    
    % Singularity warnings 
    tol = 1e-10;                 % Circular orbit tolerance 
    
    elements(5, abs(e) < tol) = zeros(1, sum(abs(e) < tol));
    elements(3, elements(4,:) == 0) = zeros(1, sum(elements(4,:) == 0));

    % Compute the semilatus rectum
    p = zeros(1,size(elements,2));
    p(1, elements(1,:) == Inf) = elements(end, elements(1,:) == Inf);                                                       % Semilatus rectum of the orbit
    p(1, elements(1,:) ~= Inf) = elements(1, elements(1,:) ~= Inf) .* (1 - elements(2,elements(1,:) ~= Inf) .^2 );          % Semilatus rectum of the orbit

    % Compute the angular momentum norm
    h = sqrt(mu .* p);                                                    % Angular momentum of the orbit
    
    % Compute the mean anomaly
    theta = Astrodynamics.KeplerSolver(e, elements(6,:));                         % True anomaly in the orbit
    
    % Compute the perifocal state vector
    r_norm = p ./ (1 + e .* cos(theta) );
    r = r_norm .* [cos(theta); sin(theta); zeros(1,size(theta,2))];       % Position vector in the perifocal frame
    v = mu ./ h .* [-sin(theta); e+cos(theta); zeros(1,size(theta,2))];   % Velocity vector in the perifocal frame

    s = zeros(6, size(theta,2));

    for i = 1:size(theta,2)
        % Rotation matrix from the inertial to the perifocal frame
        Q = Astrodynamics.euler_matrix(elements(:,i));
           
        % Output
        s(1:3,i) = Q.' * r(:,i);      % Position vector in the inertial frame
        s(4:6,i) = Q.' * v(:,i);      % Velocity vector in the inertial frame
    end
end

% Function to compute the orbital elements from the perifocal state vector
function [elements] = prv2coe(mu, s)
    % State variables 
    r = s(1:3).';                                   % Position vector
    v = s(4:6).';                                   % Velocity vector 
    
    % Compute the eccentricy and angular momentum vectors 
    h = cross(r,v);                                 % Angular momentum vector
    h = norm(h);                                    % Angular momentum norm
    
    % Compute the true anomaly 
    theta = atan2(r(2), r(1)); 
    
    % Compute the eccentricity 
    e = v(2)*(h/mu)-cos(theta); 
    
    % Compute orbital energy 
    H = norm(v)^2/2-mu/norm(r);
    
    % Determine type of orbit 
    if (e ~= 1)
        a = -mu/(2*H);                              % Semimajor axis of the orbit
        p = a*(1-e^2);                              % Semilatus rectum of the orbit
    else
        p = h^2/mu;                                 % Semilatus rectum of the orbit
        a = Inf;                                    % Semimajor axis of the orbit
    end
        
    % Compute the rest of the orbital elements
    RAAN = NaN;                                         % RAAN
    omega = NaN;                                        % Argument of perigee
    i = NaN;                                            % Inclination of the orbit
    
    % Mean anomaly
    sinE = sqrt(1-e^2)*sin(theta)/(1+e*cos(theta));     % Sine of the eccentric anomaly
    cosE = (e+cos(theta))/(1+e*cos(theta));             % Cosine of the eccentric anomaly
    E = atan2(sinE, cosE);                              % Eccentric anomaly
    M = E-norm(e)*sin(E);                               % Mean anomaly
    
    % Save the classical orbital elements 
    elements = [a e RAAN i omega M p];
    
    % Singularity warnings 
    tol = 1e-10;               % Circular orbit tolerance    
    if (abs(e) < tol)
        warning('Orbit is circular to numerical precision');   
    end
end