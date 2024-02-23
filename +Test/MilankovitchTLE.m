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
file_name = '+Test\TEST_TLE.TLE';        % Test data file name

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

fclose(fid);                             % Close file

%% Verification and validation of the propagator
% Dimensional units 
mu = tle.rec.mu;                            % Gravitational constant used        
L = tle.rec.radiusearthkm;                  % Characteristic length
T = sqrt( L^3 / mu );                       % Characteristic time                                      

% Get the initial Milankovitch elements conditions in the TEME frame from the TLE (which is approximately equivalent to Brouwer's elements)ç
meanElements = [tle.n * 2*pi/1440 tle.ecc tle.raanDeg tle.incDeg tle.argpDeg tle.maDeg];      % Mean COE elements
meanElements(3:end) = deg2rad(meanElements(3:end));                                           % All angles in radians

% Un-kozai the mean motion
meanElements = Astrodynamics.TLE2COE(tle.rec.xke, tle.rec.j2, meanElements, true);

% Initial conditions
s0 = Astrodynamics.MKV2COE(1, meanElements.', false);                         % Initial mean Milankovitch elements                                

% TEME propagated Cartesian elements (initial osculating conditions in TEME)
rv0 = tle.getRV(0);                                                         
rv0 = reshape(rv0, [], 1);
rv0(1:3,1) = rv0(1:3,1) / L;
rv0(4:6,1) = rv0(4:6,1) / L * T;

% Propagation setup 
model = 'SGP4';
options = odeset('AbsTol', 1E-22, 'RelTol', 2.24E-14);

% Elapsed time to propagate since the generation of the TLE (initial conditions)
elapsed_epoch = linspace(0, 1 * 86400, 1E2);    % Seconds since the TLE epoch

% Preallocation for speed 
RV = zeros(6,length(elapsed_epoch));            % TEME Cartesian state vector from SGP4
RV2 = RV;                                       % TEME Cartesian state vector from MKV dynamics

for i = 1:length(elapsed_epoch)

    dt = elapsed_epoch(i) - elapsed_epoch(1);   % Propagation step in seconds
    
    switch (model)
        case 'SGP4'
            % Propagate to the given epoch by means of the SGP4 model
            dt_min = dt / 60;                   % Propagation step in minutes
            rv = tle.getRV(dt_min);             % TEME propagated Cartesian elements
            rv = reshape(rv, [], 1);
            sgp_rv(1:3,1) = rv(1:3,1) / L;      % Non-dimensional position
            sgp_rv(4:6,1) = rv(4:6,1) / L * T;  % Non-dimensional velocity
    
        case 'APSO'
            % Propagate to the given epoch by means of the osculating J2 problem model
            if (dt > 0)
                tspan = linspace(0, dt / T, 1E3);
                [~, rv] = ode45(@(t,s)Astrodynamics.APSO_dynamics(1, tle.rec.j2, 1, s, 0, 0), [0 dt / T], rv0, options);
                sgp_rv = rv(end,:).';
%                 rv(1:3,1) = rv(1:3,1) * L;
%                 rv(4:6,1) = rv(4:6,1) * L / T;
            else
                sgp_rv = rv0;
            end
        
        otherwise
            error('No valid dynamics model was selected');
    end

    % Compute the very same quantities by propagating the Milankovitcch elements
    if (dt > 0)
        % Propagate them   
        tspan = linspace(0, dt / T, 1E4);
        [~, s] = ode45(@(t,s)Astrodynamics.milankovitch_dynamics(1, tle.rec.j2, [0;0;1], t, s), tspan, s0, options);
        s = s(end,:).';
    else
        s = s0;
    end

    J = Astrodynamics.milankovitch_jacobian(tle.rec.j2, [0;0;1], s);
    
    % Save the results for the further processing
    mkv_rv = Astrodynamics.MKV2ECI(1, s, true);
    LM = Astrodynamics.Lara2ECI(s0(3,1), mkv_rv, false);
    LEO = Astrodynamics.BrouwerLaraCorrections(tle.rec.j2, tle.rec.j3, s0(3,1), LM);
    mkv_rv = Astrodynamics.Lara2ECI(s0(3,1), LEO, true);

    % Save results
    RV(:,i) = sgp_rv;
    RV2(:,i) = mkv_rv;

%     RV2(1:3,i) = eci_rv(1:3,1) * L;
%     RV2(4:6,i) = eci_rv(4:6,1) * L / T;
end
        
%% Analysis
% Orbit geometry error quantities
r_error = RV(1:3,:) - RV2(1:3,:);
r_error = sqrt(dot(r_error, r_error, 1));             % Propagated angular momentum relative error to the SGP4 reference
r_error = r_error * L;
mu_rerr = mean(r_error);                              % Mean of the position vector error to the SGP4 reference
sigma_rerr = std(r_error);                            % Standard deviation of the position vector error to the SGP4 reference

v_error = RV(4:6,:) - RV2(4:6,:);
v_error = sqrt(dot(v_error, v_error, 1));             % Propagated velocity vector error to the SGP4 reference
v_error = v_error * L / T;
mu_verr = mean(v_error);                              % Mean of the velocity vector error to the SGP4 reference
sigma_verr = std(v_error);                            % Standard deviation of the velocity vector error to the SGP4 reference

%% Display
figure
subplot(2,1,1)
hold on
plot(elapsed_epoch / 60, r_error + 3 * sigma_rerr)
yline(mu_rerr, 'k--')
scatter(elapsed_epoch / 60, r_error, 'filled')
ylabel('Position error [km]')
legend('3sigma', 'Mean', 'Error')
grid on;

subplot(2,1,2)
hold on
plot(elapsed_epoch / 60, v_error + 3 * sigma_verr)
yline(mu_verr, 'k--')
scatter(elapsed_epoch / 60, v_error, 'filled')
ylabel('Velocity error [km/s]')
legend('3sigma', 'Mean', 'Error')
grid on;
xlabel('$t$ [min]')

figure 
view(3)
grid on;
hold on 
scatter3(RV(1,:),RV(2,:),RV(3,:), 'b')
scatter3(RV2(1,:),RV2(2,:),RV2(3,:), 'r')
xlabel('X')
ylabel('Y')
zlabel('Z')
legend('SGP4', 'MKV')

figure
subplot(3,1,1)
hold on
scatter(elapsed_epoch / 60, RV(1,:), 'b*')
scatter(elapsed_epoch / 60, RV2(1,:), 'r')
ylabel('X')
legend('SGP4', 'MKV')
grid on;

subplot(3,1,2)
hold on
scatter(elapsed_epoch / 60, RV(2,:), 'b*')
scatter(elapsed_epoch / 60, RV2(2,:), 'r')
ylabel('Y')
legend('SGP4', 'MKV')
grid on;

subplot(3,1,3)
hold on
scatter(elapsed_epoch / 60, RV(3,:), 'b*')
scatter(elapsed_epoch / 60, RV2(3,:), 'r')
ylabel('Z')
legend('SGP4', 'MKV')
grid on;

xlabel('$t$ [min]')