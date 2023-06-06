%% FOSSA Systems %% 
% Constellation design usign evolutionary strategies, IAC 2022 % 
% Date: 04/08/2022

%% SSO Elements
% Compute the orbital elements of the SSO as a function of the altitude 

% Inputs: - vector LTAN, the LTAN reference of the design
%         - array coe, containing Nx3 orbital elements: semimajor axis and omega

% Outputs: - frozen orbit corrected SSO orbital elements

function [COE] = sso_elements(LTAN, coe)
    % Constants 
    mu = 3.986e14;                  % Graviational parameter of the Earth 
    J2 = 1.08263e-3;                % Second zonal harmonic of the Earth
    J3 = -2.53881e-6;               % Third zonal harmonic of the Earth
    Re = 6378.14e3;                 % Mean Earth radius
    OmegaS(1) = 360/365.242199;     % Earth's angular motion round the Sun
    OmegaS(2) = deg2rad(OmegaS(1));
    OmegaS(2) = OmegaS(2)/86400;

    % LTAN sanity check 
    if (length(LTAN) == 1)
        LTAN = repmat(LTAN,size(coe,1),1);
    end

    % Preallocation 
    COE = zeros(size(coe,1),6);

    % Computation of the orbital elements 
    for i = 1:size(coe,1)
        % Setup of the iterative loop 
        tol = 1e-12;                                    % Convergence tolerance
        error = 1;                                      % Initial error

        % Initial guess
        a = coe(i,1);                                   % Semimajor axis of the orbit
        e = 0;                                          % Orbital eccentricity
        n = sqrt(mu/a^3);                               % Orbital mean motion
        p = a*(1-e^2);                                  % Semilatus rectum of the orbit
        id = acos((-2*OmegaS(2)*p^2)/(3*J2*Re^2*n));    % Desired inclination

        % Iterative solver 
        while (error >= tol)
            e = (-1/2)*(J3/J2)*(Re/a)*sin(id);                                 % Frozen eccentricity
            p = a*(1-e^2);                                                     % Semilatus rectum of the orbit
        
            % J2 perturbed mean motion
            n = n*(1+(3/2)*J2*(Re/a)^2*sqrt(1-e^2)*(1-(3/2)*sin(id)^2));       % First order perturbation in the mean motion due to the J2
            i_p = acos((-2*OmegaS(2)*p^2)/(3*J2*Re^2*n));                      % Sun-synchronous condition
        
            % Convergence criteria
            error = abs(i_p-id);
            id = i_p;
        end

        % Compute the RAAN 
        RAAN = (360/24)*LTAN(i)-180+OmegaS(1)*coe(i,2);

        % Final COE 
        COE(i,:) = [coe(i,1) e deg2rad(RAAN) id coe(i,3) 360/size(coe,1)*(i-1)];
    end
end