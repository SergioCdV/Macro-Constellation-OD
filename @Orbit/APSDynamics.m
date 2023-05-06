%% Constellation macro-orbit determination 
% Date: 01/30/2023
% Author: Sergio Cuevas del Valle

%% Artificial Satellite Problem dynamics %% 
% Function to propagate an orbit under Keplerian dynamics and J2 perturbations

% Inputs: 

% Outputs: 

function [AuxEvolution] = APSDynamics(obj, tspan, model)
    % Branch the dynamics
    switch (model)
        case 0
            % Switch the Keplerian propagation
            AuxOrbit = obj.ChangeStateFormat('ECI');
        
            % Integration of the osculating problem
            [~, AuxEvolution] = ode45(@(t,s)APSO_dynamics(obj.mu, obj.J2, obj.Re, s, obj.PropagatedEpoch, t), tspan, AuxOrbit.ElementSet, obj.IntegrationOptions);
            AuxOrbit.StateEvolution = [tspan.' AuxEvolution]; 
            AuxEvolution = AuxOrbit.ChangeStateFormat(obj.ElementType).StateEvolution(:,2:end);

        case 1
            % Switch the Keplerian propagation
            AuxElements = obj.ChangeStateFormat('COE').ElementSet;
        
            % Integration of the mean problem
            AuxEvolution = APSM_dynamics(obj.mu, obj.J2, obj.Re, AuxElements, obj.PropagatedEpoch, tspan);

        case 2
            % Preallocation 
            AuxEvolution = zeros(length(tspan), 10);

            % Integration 
            for i = 1:length(tspan)
                AuxEvolution(i,1:6) = SGP4(tspan(i), 0, obj.ElementSet);
            end

        otherwise
    end
end

%% Auxiliary variables
% J2 osculating dynamics 
function [dS] = APSO_dynamics(mu, J2, Re, s, InitialEpoch, tspan)
   r = s(1:3);      % ECI position coordinates
   v = s(4:6);      % ECI velocity coordinates 

   % Dynamics
   R = norm(r);
   gamma = -mu*r/R^3.*[1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); ...
                       1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-1); ...
                       1-3/2*J2*(Re/R)^2*(5*(r(3)/R)^2-3)];
   dS = [v; gamma];
end

% J2 mean dynamics 
function [AuxEvolution] = APSM_dynamics(mu, J2, Re, s, InitialEpoch, tspan)
    % Constants 
    a = s(1);            % Semimajor axis 
    e = s(2);            % Eccentricity 
    i = s(4);            % Inclination
    n = sqrt(mu/a^3);    % Mean motion
    p = a*(1-e^2);       % Semilatus rectum

    % Kozai's method
    dS = [-3/2*n*J2*(Re/p)^2*cos(i); ...
           3/2*n*J2*(Re/p)^2*(2-5/2*sin(i)^2); ...
          n*(1+3/2*J2*(Re/p)^2*sqrt(1-e^2)*(1-3/2*sin(i)^2))];      
    
    % Integration 
    AuxEvolution(:, [1 2 4 7]) = repmat([a e i s(7)], length(tspan), 1);
    AuxEvolution(:, [3 5 6]) = repmat(s([3 5 6]), length(tspan), 1)+(dS.*tspan).';
end

% SGP4 
function [stateVector] = SGP4(epoch, initial_epoch, TLE)
    % Model constants 
    XKMPER = 6378.135;                      % Earth radius
    aE = 1;                                 % Planet radius
    ke = 0.743669161E-1;                    % Planet gravitational parameter
    k2 = 5.413080E-4;                       % J2 perturbation
    k4 = 0.62098875E-6;                     % J4 perturbation
    J3 = -2.53881e-6;                       % J3 perturbation
    A3 = -J3*(aE^3);                        % J3 perturbation
    q0 = 1.88027916E-9;                     % Model constant
    s = 1.01222928;                         % Model constant
    rp1 = aE+98/XKMPER;                     % Perigee altitude 1
    rp2 = aE+156/XKMPER;                    % Perigee altitude 2
    rp3 = aE+220/XKMPER;                    % Perigee altitude 3
    C = zeros(1,5);                         % Model coefficients 
    D = zeros(1,4);                         % Model coefficientes
    
    % Initial orbital elements 
    n0 = TLE(1)*(2*pi/(24*60));             % Mean orbital motion
    e0 = TLE(2);                            % Eccentricity
    RAAN0 = deg2rad(TLE(3));                % Right ascension of the ascendent node
    i0 = deg2rad(TLE(4));                   % Inclination
    omega0 = deg2rad(TLE(5));               % Argument of perigee
    M0 = deg2rad(TLE(6));                   % Mean anomaly
    sB = TLE(7);                            % B* coefficient
    dt = (epoch-initial_epoch)/60;          % Time step between TLE epoch
    
    % Main computations
    a1 = (ke/n0)^(2/3);                                     % Initial semimajor axis
    d1 = (3/2)*(k2/a1^2)*((3*cos(i0)-1)/(1-e0^2)^(3/2));
    a0 = a1*(1-(1/3)*d1-d1^2-(134/81)*d1^3);                % Semimajor axis perturbation
    d0 = (3/2)*(k2/a0^2)*((3*cos(i0)-1)/(1-e0^2)^(3/2));
    ddn0 = n0/(1+d0);                                       % Original mean orbital motion
    dda0 = a0/(1-d0);                                       % Original semimajor axis
    rp = dda0*(1-e0);                                       % Perigee altitude
    
    % Compute constants 
    if (rp >= rp1) && (rp <= rp2)
        starS = dda0*(1-e0)-s-aE;
        q = (q0+s-starS)^4;
        s = starS;
    elseif (rp < rp1)
        starS = 20/(XKMPER)-aE;
        q = (q0+s-starS)^4;
        s = starS;
    else
        q = q0;
    end
    
    theta = cos(i0);
    chi = 1/(dda0-s);
    beta0 = sqrt(1-e0^2);
    eta = dda0*e0*chi;
    C(2) = q*chi^4*ddn0*(1-eta^2)^(-7/2)*(dda0*(1+(3/2)*eta^2+(4*e0)*eta+e0*eta^3)+ ...
           (3/2)*((k2*chi)/(1-eta^2))*((-1/2)+(3/2)*theta^2)*(8+24*eta^2+3*eta^4));
    C(1) = sB*C(2);
    C(3) = (q*chi^5*A3*ddn0*aE*sin(i0))/(k2*e0);
    C(4) = 2*ddn0*q*chi^4*dda0*beta0^2*(1-eta^2)^(-7/2)*((2*eta*(1+e0*eta)+ ...
           (1/2)*e0+(1/2)*eta^3)+((-2*k2*chi)/(dda0*(1-eta^2)))* ...
           (3*(1-3*theta^2)*(1+(3/2)*eta^2-2*e0*eta-(1/2)*e0*eta^3)+ ...
           (3/4)*(1-theta^2)*(2*eta^2-e0*eta-e0*eta^3)*cos(2*omega0)));
    C(5) = 2*q*chi^4*dda0*beta0^2*(1-eta^2)^(-7/2)*(1+(11/4)*eta*(eta+e0)+e0*eta^3);
    D(2) = 4*dda0*chi*C(1)^2;
    D(3) = (4/3)*dda0*chi^2*(17*dda0+s)*C(1)^3;
    D(4) = (2/3)*dda0*chi^3*(221*dda0+31*s)*C(1)^4;
    
    % Secular effects of drag and gravitation 
    MDF = M0+(1+((3*k2*(3*theta^2-1))/(2*dda0^2*beta0^3))+ ...
                ((3*k2^2*(13-78*theta^2+137*theta^4)/(16*dda0^4*beta0^7))))*ddn0*dt;
    omegaDF = omega0+(((3*k2^2*(7-114*theta^2+395*theta^4))/(16*dda0^4*beta0^8))+ ... 
                      ((5*k4*(3-36*theta^2+49*theta^4))/(4*dda0^4*beta0^8))+ ...
                      ((-3*k2*(1-5*theta^2))/(2*dda0^2*beta0^4)))*ddn0*dt;
    RAANDF = RAAN0+(((3*k2^2*(4*theta-19*theta^3))/(2*dda0^4*beta0^8))+ ...
                    ((5*k4*theta*(3-7*theta^2))/(2*dda0^4*beta0^8))+ ...
                    ((-3*k2*theta)/(dda0^2*beta0^4)))*ddn0*dt;
    
    if (rp <= rp3)
        domega = 0; 
        dM = 0;
        C(5) = 0;
    else
        domega = sB*C(3)*cos(omega0)*dt;
        dM = (-2/3)*q*sB*chi^4*(aE/(e0*eta))*((1+eta*cos(MDF))^3-(1+eta*cos(M0))^3);
    end
    
    Mp = MDF+dM+domega;
    omega = omegaDF+domega-dM;
    RAAN = RAANDF-(21/2)*((ddn0*k2*theta)/(dda0^2*beta0^2))*C(1)*dt^2;
    e = e0-sB*C(5)*(sin(Mp)-sin(M0))-sB*C(4)*dt;
    
    if (rp <= rp3)
        a = dda0*(1-C(1)*dt)^2;
        L = Mp+omega+RAAN+ddn0*((3/2)*C(1)*dt^2);
    else
        a = dda0*(1-C(1)*dt-D(2)*dt^2-D(3)*dt^3-D(4)*dt^4)^2;
        L = Mp+omega+RAAN+ddn0*((3/2)*C(1)*dt^2+(D(2)+2*C(1)^2)*dt^3+ ... 
            (1/4)*(3*D(3)+12*C(1)*D(2)+10*C(1)^3)*dt^4+ ...
            (1/5)*(3*D(4)+12*C(1)*D(3)+6*D(2)^2+30*C(1)^2*D(2)+15*C(1)^4)*dt^5);
    end
    
    beta = sqrt(1-e^2);     % Semilatus rectum parameter
    n = ke/a^(3/2);         % Mean orbital motion
       
    % Long periodic terms 
    axN = e*cos(omega);
    Ll = ((A3*sin(i0))/(8*k2*a*beta^2))*(e*cos(omega))*((3+5*theta)/(1+theta));
    ayNL = (A3*sin(i0))/(4*k2*a*beta^2);
    ayN = e*sin(omega)+ayNL;
    Lt = Ll+L;
    trueE = kepler_spg4(Lt, RAAN, axN, ayN);
    ecE = axN*cos(trueE)+ayN*sin(trueE);
    esE = axN*sin(trueE)-ayN*cos(trueE);
    eL = sqrt(axN^2+ayN^2);
    pL = a*(1-eL^2);
    r = a*(1-ecE);
    dr = ke*(sqrt(a)/r)*esE;
    rdf = ke*sqrt(pL)/r;
    cu = (a/r)*(cos(trueE)-axN+(ayN*esE)/(1+sqrt(1-eL^2)));
    su = (a/r)*(sin(trueE)-ayN-(axN*esE)/(1+sqrt(1-eL^2)));
    u = atan2(su, cu); 
    deltaR = (k2/(2*pL))*(1-theta^2)*cos(2*u);
    deltaU = -(k2/(2*pL)^2)*(7*theta^2-1)*sin(2*u);
    deltaRAAN = (k2/(2*pL))*(theta)*sin(2*u);
    deltaI = (3*k2/(2*pL^2))*(sin(i0))*sin(2*u);
    deltadr = (-k2*n/pL)*(1-theta^2)*sin(2*u);
    deltardf = (k2*n/pL)*((1-theta^2)*cos(2*u)-3/2*(1-3*theta^2));
    r = deltaR+r*(1-(3/2)*k2*((sqrt(1-eL^2))/pL^2)*(3*theta^2-1));
    u = u+deltaU;
    RAAN = RAAN+deltaRAAN; 
    i = i0+deltaI; 
    dr = dr+deltadr;
    rdf = rdf+deltardf;
    M = [-sin(RAAN)*cos(i); cos(omega)*cos(i); sin(i)];
    N = [cos(RAAN); sin(RAAN); 0];
    U = M*sin(u)+N*cos(u);
    V = M*cos(u)-N*sin(u);
     
    % Output 
    r = r*U;                    % Position vector
    r = XKMPER*1e3*r;           % Dimensinal position vector
    v = dr*U+rdf*V;             % Velocity vector
    v = (XKMPER*1e3/60)*v;      % Dimensional velocity vector
    stateVector = [r; v];       % Final state vector 
end

%% Auxiliary functions 
% Solve Kepler equation for the true eccentric longitude for SGP4
function [trueE] = kepler_spg4(Lt, RAAN, axN, ayN)
    %Constants 
    tol = 1e-3;          % Newton method tolerance
    maxIter = 1e4;      % Maximum number of iterations
    iter = 1;           % Initial iteration
    go_on = true;       % Convergence flag
    U = Lt-RAAN;        % Initial true eccentric longitude
    E(iter) = U;        % Initial true eccentric longitude to iterate    
    
    % Main loop 
    while ((go_on) && (iter < maxIter))
        % Main computations 
        dE = (U-ayN*cos(E(iter))+axN*sin(E(iter))-E(iter))/(1-axN*cos(E(iter))-ayN*sin(E(iter)));
        E(iter+1) = E(iter)+dE;
        
        % Convergence checking
        if (abs(E(iter+1)-E(iter)) <= tol)
            go_on = false;
        else
            iter = iter+1;
        end
    end

    % Output 
    trueE = E(end);     % Final true eccentric longitude
end

