%% Orbital estimation Example %%
% Date: 22/01/2024
% This file provides an UKF filter implementation for the Milankovitch dynamics

clear; 
close all; 
clc;
rng(1); 
format long; 

set_graphics();

%% Initial conditions and set up 
% Canonical units are assumed for the sake of simplicity
Re = 1;                 % Radius of the Earth (in Earth radii)
J2 = 1.08263e-3;        % J2 parameter of the Earth
J3 = J2^2;                  % J3 parameter of the Earth
Keci = [0;0;1];         % Third unit vector in the ECI reference frame triad

% Milankovitch initial conditions in the ECI frame (non-dimensional units)
h0 = [0; 0; 1.8];       % Initial angular momentum vector (equatorial orbit) [-]
e0 = [0.1; 0; 0];     % Initial eccentricity vector (equatorial orbit of 0.001 eccentricity, quasi-circular) [-]
l0 = 0;                 % Initial longitude [rad]
s0 = [h0; e0; l0];      % Complete initial state

% Integration setup 
options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14);      % Integration tolerances

% Dynamics time span 
N = 1E3;                        % Number of steps
t0 = 0;                         % Initial epoch
tf = 20;                        % Final epoch [orbital periods at the Earth's WGS84 radius]
tspan = linspace(t0, tf, N);    % Integration time span


%% Propagation 
[tp, s] = ode45(@(t,s)utils.milankovitch_dynamics(Re, J2, Keci, t, s), tspan, s0, options);

% Reduced MKV state 
s = s(:,[1:5 7]);       % Disregard the inertial z component of the eccentricity vector

%% Lara's transformation test 
% Determine the third component of the eccentricity to handle a 6 x 1 state vector
e = -dot(s(1,1:2), s(1,4:5), 1) ./ s(1,3);
st = [s(1,1:5) e s(1,6)].';

% Transform to osculating position and velocity 
mkv_rv = utils.MKV2ECI(1, st, true);                                        % From MKV to Cartesian
LM = utils.Lara2ECI(st(3,1), mkv_rv, false);                                % From Cartesian to Lara
LEO = utils.BrouwerLaraCorrections(1, J2, J3, st(3,1), LM);                 % Mean to oculating Lara
rv = utils.Lara2ECI(st(3,1), LEO, true);                                    % Osculating Cartesian
LEOb = utils.Lara2ECI(st(3,1), rv, false);                                  % Osculating Lara

test(1) = norm(LEO-LEOb);

LMb = utils.BrouwerLaraCorrections(1, J2, J3, st(3,1), LEOb, 0);            % Mean Lara
test(2) = norm(LM-LMb);

%% Offline observation 
% Constants 
Pobs = blkdiag( 1e-9 * eye(3), 1e-9 * eye(3), 1E-2 );             % Covariance of the observation 

% Assume a linear identity observation model (H = I)
PD = 0.3;                                                        % Probability of observation
index = logical(randsrc(size(s,1), 1, [0, 1; 1-PD, PD]));        % Observation epochs

% Observed state vector
sobs = [tp(index) s(index,:)];

meas = zeros(size(sobs,1), 8);                              % Preallocation of the measurements
meas(:,1) = tp(index);
for i = 1:size(sobs,1)
    meas(i,2:end) = observation_model(sobs(i,2:end).').';   % Go through the measurement model
    meas(i,2:end) = meas(i,2:end) + mvnrnd(zeros(size(meas(i,2:end))), Pobs, 1);         % Add some noises
end

%% Offline estimation (EKF)
fprintf('Initialization... \n');
fprintf('------------------------------\n');
fprintf('Running the filter: \n');

% Constants 
Q = 1e-5 * eye(6);      % Dynamical noise 
R = 1e-3 * eye(7);      % Measurement noise 

s_dim = 6;              % Dimension of the (reduced) state vector
y_dim = 7;              % Dimension of the observation vector

% Estimator setup
beta = 2;               % Spread of the sigma points
alpha = 1E-3;           % Spread of the sigma points
k = 0;                  % Spread of the sigma points

UKF_estimator = aux_filters.UKF('UKF-S', beta, alpha, k);
UKF_estimator = UKF_estimator.AssignStateProcess(s_dim, @(s, time_step)state_model(1, J2, Keci, s, time_step));
UKF_estimator = UKF_estimator.AdditiveCovariances(Q, R);
UKF_estimator = UKF_estimator.AssignObservationProcess(y_dim, @(s)observation_model(s));
UKF_estimator = UKF_estimator.Init();

% Initial conditions 
se = [0; 0; 1.5; 0; 0; 0];        % Estimator initial conditions 
Pe = 1e-3 * eye(6);               % Initial covariance matrix

% Main loop
i = 1;
epoch = tp(1);                    % Initial epoch

while (i <= size(meas,1))
    % Propagation 
    tic 
    dt = meas(i,1) - epoch;
    UKF_estimator.State = se;
    UKF_estimator.Sigma = Pe;
    [sigma_points, se, Pe] = UKF_estimator.PropagationStep(dt);

    % Correction
    [se, Pe, ~] = UKF_estimator.CorrectionStep(sigma_points, se, Pe, meas(i,2:end).');

    % Save results
    S(:,i) = se;                                   % Corrected state
    P(:,i) = reshape(Pe, [], 1);                   % Corrected covariance
    E(i) = 0.5 * log(det(2*pi * exp(1) * Pe));     % Entropy characterization 

    time(i) = toc;
    fprintf('Iteration running time: %.4f s\n', time(i));

    % Update the timer
    epoch = meas(i,1);
    i = i + 1;
end

fprintf('------------------------------\n');
fprintf('EKF filter recursion finished.\n');
fprintf('Total running time: %.4f s\n', sum(time));

%% Error estimation 
% Final observation time
t = meas(:,1);

% Compute the errors 
e = sobs(:,2:end).' - S;

% Compute the 3sigma bounds and rotate the output reference frame to diagonalize P
sigma = zeros(s_dim, length(t));
% for i = 1:length(t)
%     Sigma = reshape(P(:,i), s_dim, s_dim);
%     [V, lambda] = eig(Sigma);
%     
%     sigma(:,i) = real( sqrt( diag(lambda) ) );
%     e(:,i) = V \ e(:,i);
% end

error = e.';
sigma = sigma.';

%% Plot and results
% Propagation of the state vector and measurements 
figure 
subplot(3,1,1)
hold on
plot(tp, s(:,1:3))
scatter(meas(:,1), meas(:,2), 'Marker', 'x')
scatter(meas(:,1), meas(:,3), 'Marker', 'x')
scatter(meas(:,1), meas(:,4), 'Marker', 'x')
hold off
grid on
ylabel('$\|\mathbf{h}\|$ [-]')
legend('$h_x$', '$h_y$', '$h_z$')

e = -dot(s(:,1:2), s(:,4:5), 2) ./ s(:,3);
subplot(3,1,2)
hold on
plot(tp, [s(:,4:5) e])
scatter(meas(:,1), meas(:,5), 'Marker', 'x')
scatter(meas(:,1), meas(:,6), 'Marker', 'x')
scatter(meas(:,1), meas(:,7), 'Marker', 'x')
hold off
grid on
ylabel('$\|\mathbf{e}\|$ [-]')
legend('$e_x$', '$e_y$', '$e_z$')

subplot(3,1,3)
hold on
plot(tp, s(:,6), 'b')
scatter(meas(:,1), meas(:,8), 'r', 'Marker', 'x')
hold off
ylabel('$l$ [rad]')
grid on
xlabel('$t$')

% Conservation quantities (V&V purposes)
figure
hold on
subplot(1,2,1)
plot(tp, s(:,3), 'b')
ylabel('$h_z$ [-]')
grid on
subplot(1,2,2)
plot(tp, sqrt( dot([s(:,4:5) e], [s(:,4:5) e], 2) ), 'b')
ylabel('$\|\mathbf{e}\|$ [-]')
grid on
xlabel('$t$')

% Error plots 
figure 
subplot(3,1,1)
hold on
plot(t, error(:,1), 'b')
plot(t, error(:,1) + 3 * sigma(:,1))
plot(t, error(:,1) - 3 * sigma(:,1))
grid on
ylabel('$\Delta h_x$ [-]')
subplot(3,1,2)
hold on
plot(t, error(:,2), 'b')
plot(t, error(:,2) + 3 * sigma(:,2))
plot(t, error(:,2) - 3 * sigma(:,2))
grid on
ylabel('$\Delta h_y$ [-]')
subplot(3,1,3)
hold on
plot(t, error(:,3), 'b')
plot(t, error(:,3) + 3 * sigma(:,1))
plot(t, error(:,3) - 3 * sigma(:,1))
grid on
ylabel('$\Delta h_z$ [-]')
xlabel('$t$')

figure 
subplot(2,1,1)
plot(t, error(:,4), 'b')
plot(t, error(:,4) + 3 * sigma(:,4))
plot(t, error(:,4) - 3 * sigma(:,4))
grid on
ylabel('$\Delta e_x$ [-]')
subplot(2,1,2)
plot(t, error(:,5), 'b')
plot(t, error(:,5) + 3 * sigma(:,5))
plot(t, error(:,5) - 3 * sigma(:,5))
grid on
ylabel('$\Delta e_y$ [-]')
% subplot(3,1,3)
% plot(t, error(:,6), 'b')
% plot(t, error(:,6) + 3 * sigma(:,6))
% plot(t, error(:,6) - 3 * sigma(:,6))
% grid on
% ylabel('$\Delta e_z$ [-]')
% xlabel('$t$')

figure 
hold on
plot(t, error(:,6), 'b')
plot(t, error(:,6) + 3 * sigma(:,6))
plot(t, error(:,6) - 3 * sigma(:,6))
grid on
ylabel('$\Delta l$ [rad]')
xlabel('$t$')

figure 
hold on
plot(t, E, 'b')
grid on
ylabel('det $P$ [-]')
xlabel('$t$')

%% Auxiliary functions 
% Propagation model 
function [s] = state_model(Re, J2, Keci, s0, dt)
    % Integration of the state vector 
    options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14);      % Integration tolerances

    % Determine the third component of the eccentricity to handle a 6 x 1 state vector
    e = -dot(s0(1:2,:), s0(4:5,:), 1) ./ s0(3,:);
    s0 = [s0(1:5,:); e; s0(6,:)];

    % Integrate the whole model
    if (dt > 0)
        S0 = reshape(s0, [], 1);
        [~, saux] = ode45(@(t,s)utils.milankovitch_dynamics(Re, J2, Keci, t, s), [0 dt], S0, options);
        s = saux(end,:); 
        s = reshape(s, size(s0));
    else
        s = s0;
    end

    % Eliminate the third component of the eccentricity to handle a 6 x 1 state vector
    s = s([1:5 7],:);
end

% Observation model 
function [y, H] = observation_model(s)
    % Determine the third component of the eccentricity to handle a 6 x 1 state vector
    e = -dot(s(1:2,:), s(4:5,:), 1) ./ s(3,:);
    s = [s(1:5,:); e; s(6,:)];
    
    % Linear identity model
    y = s;                                   % Measurements
    H = eye(size(s,1), size(s,1));           % Observation matrix
end

% RADAR observation model
function [y, H] = radar_model(structData, PARAMS, init_epoch, t, s)
    % Determine the third component of the eccentricity to handle a 6 x 1 state vector
    e = -dot(s(1:2,:), s(4:5,:), 1) ./ s(3,:);
    s = [s(1:5,:); e; s(6,:)];

    % Transform to osculating position and velocity 
    mkv_rv = utils.MKV2ECI(1, s, true);                                                 % From MKV to Cartesian
    LM = utils.Lara2ECI(s(3,1), mkv_rv, false);                                         % From Cartesian to Lara
    LEO = utils.BrouwerLaraCorrections(1, tle.rec.j2, tle.rec.j3, s0(3,1), LM);         % Mean to oculating Lara
    rv = utils.Lara2ECI(s(3,1), LEO, true);                                             % Osculating Cartesian  

    % Apply the radar model 
    dateObs = init_epoch + t;
    y = observationRadar(structData.xStateatEpoch', dateObs, rv(1:3,1)*1e-3, PARAMS);  % Measurements

    H = [];           % Observation matrix (not neede for UKF)
end

function set_graphics()
    %Set graphical properties
   set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
    set(groot, 'defaultAxesFontSize', 11); 
    set(groot, 'defaultAxesGridAlpha', 0.3); 
    set(groot, 'defaultAxesLineWidth', 0.75);
    set(groot, 'defaultAxesXMinorTick', 'on');
    set(groot, 'defaultAxesYMinorTick', 'on');
    set(groot, 'defaultFigureRenderer', 'painters');
    set(groot, 'defaultLegendBox', 'off');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultLegendLocation', 'best');
    set(groot, 'defaultLineLineWidth', 1); 
    set(groot, 'defaultLineMarkerSize', 3);
    set(groot, 'defaultTextInterpreter','latex');
end
