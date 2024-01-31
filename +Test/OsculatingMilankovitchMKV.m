%% Milankovitch vs SGP4 %%
% Author: Sergio Cuevas 
% Date: 12/12/2023

% This script aims to verify the performance of the Milankovitch elements
% propagator under J2 perturbation against the standard SGP4 model. 
% 
% To do so, the osculating, propagated TLE data are transformed to the
% corresponding angular momentum and eccentricity vectors, realised in the
% TEME frame. These are then compared to those obtained by means of the
% Milankovitch propagator. 
% 
% Error performance is then analysed statistically in terms of the relative
% error between the two models obtained for the angular momentum magnitude,
% the eccentricity vector norm and the attitude error between the perifocal
% reference frames. Note that if the two models are exactly equivalent, all
% erros should be zero.

clear; 
clc;

%% Input data
file_name = '+Test\SGP4-VER.TLE';        % Test data file name
output_name = '+Test\tcppver.out';       % Output verification data

fid = fopen(file_name);                  % Open file
tlecnt = 1;                              % TLE counter
tles = cell(tlecnt);                     % Preallocation of a TLE cell array structure
line1 = fgetl(fid);

while ( ischar(line1) )
    if (line1(1) == '1')
        line2 = fgetl(fid);
        tle = SGP4.TLE(line1, line2);    % Create the TLE object
        tles{tlecnt} = tle;              % Save the TLE object for its later use
        tlecnt = tlecnt + 1;
    end
    line1 = fgetl(fid);
end

fclose(fid);                    % Close file

model = 'SGP4';
options = odeset('AbsTol', 1E-22, 'RelTol', 2.24E-14);

%% Verification and validation of the propagator
% Preallocation of variables
curobj = 0;     % Current object index
strind = 0;
prevrv = 0;
ind = 0;
cnt = 0;

% Loop over the output verification file to compute Cartesian elements for various TLEs and epochs
fid = fopen(output_name);
line1 = fgetl(fid);

while ( ischar(line1) )
    strind = strfind(line1, 'xx'); 

    if (strind > 0)
        % Get the object number
        curobj = str2double(strtrim(line1(1:strind-1)));

        % Find tle with the matching object number
        ind = 1;
        while (ind <= tlecnt && tles{ind}.objectNum ~= curobj)
            ind = ind+1;
        end
        tle = tles{ind};

    elseif (curobj > 0)
        %% PROPAGATION
        % Initial conditions for the Milankovitch elements
        mu = tle.rec.mu;                            % Gravitational constant used (check for correct units in the other propagator!!!)
        L = tle.rec.radiusearthkm;                  % Characteristic length
        T = sqrt( L^3 / mu );                       % Characteristic time 

        % Get the initial Milankovitch elements conditions in the TEME frame from the TLE (which is approximately equivalent to Brouwer's elements)
        a = (mu / (2*pi/86400 * tle.n)^2)^(1/3);                                    % Mean semimajor axis
        meanElements = [a tle.ecc tle.raanDeg tle.incDeg tle.argpDeg tle.maDeg];    % Mean COE elements
        meanElements(3:end) = deg2rad(meanElements(3:end));                         % All angles in radians
        meanElements(1) = meanElements(1) / L;                                      % Non-dimensional semimajor axis
        s0 = Astrodynamics.ECI2COE(1, meanElements, false).';                       % Cartesian elements
        s0 = Astrodynamics.MKV2ECI(1, s0, false);                                   % Initial Milankovitch elements                                  

        % Propagation
        verout = sscanf(line1, '%f');   % TLE output data for verification purposes
        elapsed_epoch = verout(1);      % Elapsed time to propagate since the generation of the TLE (initial conditions)

        % Propagate to the given epoch by means of the SGP4 model
        switch (model)
            case 'SGP4'
                % Propagate to the given epoch by means of the SGP4 model
                rv = tle.getRV(elapsed_epoch);                  % TEME propagated Cartesian elements
                rv = reshape(rv, [], 1);

            case 'APSO'
                % Propagate to the given epoch by means of the osculating J2 problem model
                rv0 = tle.getRV(0);                             % TEME propagated Cartesian elements (initial osculating conditions)
                rv0 = reshape(rv0, [], 1);
                
                % Propagation
                if (elapsed_epoch > 0)
                    tspan = linspace(0, elapsed_epoch * 60, 1E1) / T;
                    rv0(1:3,1) = rv0(1:3,1) / L;
                    rv0(4:6,1) = rv0(4:6,1) / L * T;
                    [~, rv] = ode113(@(t,s)Astrodynamics.APSO_dynamics(1, tle.rec.j2, 1, s, 0, 0), tspan, rv0);
                    rv = rv(end,:).';
                    rv(1:3,1) = rv(1:3,1) * L;
                    rv(4:6,1) = rv(4:6,1) * L / T;
                else
                    rv = rv0;
                end
            
            otherwise
                error('No valid dynamics model was selected');
        end

        % Original code reused RV so if there is an error we need to carry old values
        if (tle.sgp4Error > 0)
            rv = prevrv;
        end

        % Compute the osculating, propagated SPG4 perifocal triad (angular momentum, eccentricity and Hamilton's vector) in the TEME frame
        cnt = cnt + 1;
        h(:,cnt) = cross(rv(1:3,1), rv(4:6,1));                                     % Specific angular momentum
        e(:,cnt) = cross(rv(4:6,1), h(:,cnt)) / mu - rv(1:3,1) / norm(rv(1:3,1));   % Eccentricity vector
        y(:,cnt) = cross(h(:,cnt), e(:,cnt));                                       % Hamilton's vector (semiminor axis vector)

        %% Compute the very same quantities by propagating the Milankovitcch elements
        if (elapsed_epoch > 0)
            % Propagate them 
            elapsed_epoch = elapsed_epoch * 60;             % Propagation time in seconds   
            tspan = linspace(0, elapsed_epoch, 1E2) / T;
            [~, s] = ode113(@(t,s)Astrodynamics.milankovitch_dynamics(tle.rec.j2, [0;0;1], t, s), tspan, s0, options);
            s = s(end,:).';
        else
            s = s0;
        end

        % Save the results for the further processing
        hm(:,cnt) = s(1:3,1) * (L^2/T);             % This shall be replaced by the Milankovitch propagation and osculating transformation
        em(:,cnt) = s(4:6,1);                       % This shall be replaced by the Milankovitch propagation and osculating transformation
        ym(:,cnt) = cross(hm(:,cnt), em(:,cnt));    % This shall be replaced by the Milankovitch propagation and osculating transformation

        prevrv = rv;
    end

    line1 = fgetl(fid);
end
fclose(fid);

%% Analysis
% Orbit geometry error quantities
h_error = h - hm;
h_error = sqrt(dot(h_error, h_error, 1));             % Propagated angular momentum relative error to the SGP4 reference
mu_herr = mean(h_error);                              % Mean of the angular momentum relative error to the SGP4 reference
sigma_herr = std(h_error);                            % Standard deviation of the angular momentum relative error to the SGP4 reference

e_error = e - em;
e_error = sqrt(dot(e_error, e_error, 1));             % Propagated angular momentum relative error to the SGP4 reference
mu_eerr = mean(e_error);                              % Mean of the angular momentum relative error to the SGP4 reference
sigma_eerr = std(e_error);                            % Standard deviation of the angular momentum relative error to the SGP4 reference

% Perifocal error quantities 
q_error = zeros(1,length(h_error));

for i = 1:cnt 
    P = [e(:,i) / norm(e(:,i)) y(:,i) / norm(y(:,i)) h(:,i) / norm(h(:,i))]; 
    Pm = [em(:,i) / norm(em(:,i)) ym(:,i) / norm(ym(:,i)) hm(:,i) / norm(hm(:,i))];
    dP = P * Pm.';
    cos_theta = 0.5 * (trace(dP)-1);
    q_error(i) = rad2deg( acos( sign(cos_theta) * max(abs(cos_theta), 1) ) );
end

q_error = unwrap(q_error); 

%% Display
figure
subplot(2,1,1)
hold on
plot(1:cnt, h_error + 3 * sigma_herr)
yline(mu_herr, 'k--')
scatter(1:cnt, h_error, 'filled')
ylabel('Angular momemtum error [%]')
legend('3sigma', 'Mean', 'Error')
grid on;

subplot(2,1,2)
hold on
plot(1:cnt, e_error + 3 * sigma_eerr)
yline(mu_eerr, 'k--')
scatter(1:cnt, e_error, 'filled')
ylabel('Eccentricity error [%]')
legend('3sigma', 'Mean', 'Error')
grid on;
xlabel('TLE case')

figure
hold on
plot(1:cnt, q_error + 3 * std(q_error))
yline(mean(q_error), 'k--')
scatter(1:cnt, q_error, 'filled')
ylabel('Perifocal attitude error [deg]')
legend('3sigma', 'Mean', 'Error')
grid on;
xlabel('TLE case')