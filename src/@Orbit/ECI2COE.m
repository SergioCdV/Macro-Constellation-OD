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

function [elements] = ECI2COE(obj, s, direction)
    % Constants
    mu = obj.mu; 
    
    % Main computation 
    if (direction)
        % Preallocation 
        elements = zeros(size(s,1),7); 

        for i = 1:size(s,1)
            elements(i,:) = rv2coe(mu,s(i,:));
        end
    else
        % Preallocation 
        elements = zeros(size(s,1),6); 

        for i = 1:size(s,1)
            elements(i,:) = coe2state(mu,s(i,:));
        end
    end
end

%% Auxiliary functions 
% Function to compute the orbital elements from the inertial state vector 
function [elements] = rv2coe(mu, s)
    % State variables 
    s = s.';
    
    r = s(1:3);                                     % Position vector
    v = s(4:6);                                     % Velocity vector 
    
    % Compute the eccentricy and angular momentum vectors 
    h = cross(r,v);                                 % Angular momentum vector
    e = cross(v,h)/mu-r/norm(r);                    % Eccentricity vector
    K = [0; 0; 1];                                  % Inertial Z axis unit vector
    n = cross(K, h);                                % Node vector
    
    % Compute orbital energy 
    H = norm(v)^2/2-mu/norm(r);
    
    % Determine type of orbit 
    if (norm(e) ~= 1)
        a = -mu/(2*H);                              % Semimajor axis of the orbit
        p = a*(1-norm(e)^2);                        % Semilatus rectum of the orbit
    else
        p = norm(h)^2/mu;                           % Semilatus rectum of the orbit
        a = Inf;                                    % Semimajor axis of the orbit
    end
    
    % Compute the unit perifocal triad 
    i = e/norm(e); 
    k = h/norm(h); 
    j = cross(k,i);
    
    % Compute the rotation matrix
    Q = [i.'; j.'; k.'];                            % Perifocal rotation matrix
    r0 = Q*r;                                       % Position in the perifocal frame 
    
    % Compute the rest of the orbital elements
    RAAN = atan2(Q(3,1),-Q(3,2));                                        % RAAN
    omega = atan2(Q(1,3),Q(2,3));                                        % Argument of perigee
    i = acos(Q(3,3));                                                    % Inclination of the orbit
    
    % Mean anomaly
    theta = atan2(r0(2), r0(1));                                         % True anomaly of the orbit
    sinE = sqrt(1-norm(e)^2)*sin(theta)/(1+norm(e)*cos(theta));          % Sine of the eccentric anomaly
    cosE = (norm(e)+cos(theta))/(1+norm(e)*cos(theta));                  % Cosine of the eccentric anomaly
    E = atan2(sinE, cosE);                                               % Eccentric anomaly
    M = E-norm(e)*sin(E);                                                % Mean anomaly
        
    % Save the classical orbital elements 
    elements = [a norm(e) RAAN i omega M p];
    elements = rv_singularity(e, n, r, Q, elements);
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

% Function to handle orbital elements singularities when converting from the inertial state vector
function [elements] = rv_singularity(ev, n, r, Q, elements)
    % State variables 
    e = elements(2);           % Orbit eccentricity
    i = elements(4);           % Orbit inclination 
    M = elements(6);           % Mean anomaly
    tol = 1e-10;               % Circular orbit tolerance
    
    % Singularity warnings 
    if (any(isnan(Q)))
        warning('Euler angles are numerically ill-conditioned');
    end
    
    if (abs(e) < tol)
        warning('Orbit is circular to numerical precision');
    end
    
    if (abs(i) < tol)
        warning('Orbit is equatorial to numerical precision');
    end
    
    %Singularity handling 
    if (abs(e) < tol)
        if (abs(i) < tol)
            M = acos(r(1)/norm(r));                 % Circular equatorial orbit
            if (r(2) < 0)
                M = 2*pi-M;
            end
        else
            M = acos(dot(n,r)/(norm(n)*norm(r)));   % Circular inclined orbit
            if (r(3) < 0)
                M = 2*pi-M;
            end
        end
    else
        if (abs(i) < tol)
            M = acos(ev(1)/norm(ev));               % Equatorial orbit
            if (ev(2) < 0)
                M = 2*pi-M; 
            end
        end
    end
    
    % Reconstruction of the orbital elements 
    elements(6) = M;
end

% Transform COE to Cartesian state vector
% Inputs: - scalar mu, the gravitational parameter of the central voyd.
%         - vector elements, containing the classical orbital elements. 
% Ouputs: - vector s, containing the state vector in the inertial frame (position and velocity vectors).
function [s] = coe2state(mu, elements)
    % Constants 
    e = elements(2);                                            % Eccentricity of the orbit
    
    % Singularity warnings 
    tol = 1e-10;                                                % Circular orbit tolerance    
    if (abs(norm(e)) < tol)
        elements(5) = 0;
    end
    
    if (elements(4) == 0)
        elements(3) = 0;
    end
    
    % Compute the semilatus rectum
    if (elements(1) == Inf) 
        p = elements(end);                                      % Semilatus rectum of the orbit
    else
        p = elements(1)*(1-elements(2)^2);                      % Semilatus rectum of the orbit
    end
    
    % Compute the angular momentum norm
    h = sqrt(mu*p);                                             % Angular momentum of the orbit
    
    % Compute the mean anomaly
    theta = kepler(elements);                                   % True anomaly in the orbit
    
    % Compute the perifocal state vector
    r = (p/(1+e*cos(theta)))*[cos(theta); sin(theta); 0];       % Position vector in the perifocal frame
    v = mu/h*[-sin(theta); e+cos(theta); 0];                    % Velocity vector in the perifocal frame
    
    % Rotation matrix from the inertial to the perifocal frame
    Q = euler_matrix(elements);
       
    % Output
    r = Q.'*r;      % Position vector in the inertial frame
    v = Q.'*v;      % Velocity vector in the inertial frame
    s = [r; v];     % State vector
end

% Kepler equation 
function [theta] = kepler(elements)
    % Constants 
    e = elements(2);                                            % Eccentricity of the orbit
    M0 = elements(6);                                           % Mean anomaly of the orbit 
    
    % Set up the loop 
    k = 5;                                                      % Conway constant
    tol = 1e-15;                                                % Convergence tolerance
    iterMax = 10^6;                                             % Maximum number of iterations
    GoOn = true;                                                % Convergence flag
    iter = 1;                                                   % Initial iteration
    u = M0+e;                                                   % Conway method variable
    E(iter) = (M0*(1-sin(u))+u*sin(M0))/(1+sin(M0)-sin(u));     % Initial guess for the eccentric anomaly
    
    % Main computation 
    while ((GoOn) && (iter < iterMax))
        % Laguerre-Conway iterations
        f = E(iter)-e*sin(E(iter))-M0; 
        df = 1-e*cos(E(iter));
        ddf = e*sin(E(iter));
        dn = -k*(f)/(df+sqrt((k-1)^2*df^2-k*(k-1)*f*ddf));
        E(iter+1) = E(iter)+dn;
        
        % Convergence checking 
        if (abs(dn) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end  
    
    % True anomaly
    theta = atan2(sqrt(1-e^2)*sin(E(end))/(1-e*cos(E(end))), (cos(E(end))-e)/(1-e*cos(E(end))));
end

% ECI to PF rotation matrix
% Inputs: - vector elements, the mean classical Euler elements
% Outputs: - matrix Q, the rotation matrix from the inertial to the perifocal frame
function [Q] = euler_matrix(elements)
    % Elements of interest 
    RAAN = elements(3); 
    i = elements(4); 
    omega = elements(5); 
    
    % Compute the rotation matrix (Euler sequence ZXZ)
    Q1 = [cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];
    Q2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    Q3 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
    Q = Q3*Q2*Q1;
end