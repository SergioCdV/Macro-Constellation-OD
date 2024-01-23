%% Orbital estimation Example %%
% Author: Sergio Cuevas 
% Date: 22/01/2024
% This file provides an EKF filter implementation for the Milankovitch dynamics

clear; 
close all; 
clc;
rng(1); 
format long; 
set_graphics(); 

%% Initial conditions and set up 
% Canonical units are assumed for the sake of simplicity
J2 = 1.08263e-3;        % J2 parameter of the Earth
Keci = [0;0;1];         % Third unit vector in the ECI reference frame triad

% Milankovitch initial conditions in the ECI frame
h0 = [0; 0; 1.8];         % Initial angular momentum vector (equatorial orbit) [-]
e0 = [0.001; 0; 0];     % Initial eccentricity vector (equatorial orbit of 0.001 eccentricity, quasi-circular) [-]
l0 = 0;                 % Initial longitude [rad]
s0 = [h0; e0; l0];      % Complete initial state

% Integration setup 
options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14);      % Integration tolerances

% Dynamics time span 
N = 1E3;                        % Number of steps
t0 = 0;                         % Initial epoch
tf = 20;                         % Final epoch [orbital periods at the Earth's WGS84 radius]
tspan = linspace(t0, tf, N);    % Integration time span

%% Propagation 
[tp, s] = ode45(@(t,s)Astrodynamics.milankovitch_dynamics(J2, Keci, t, s), tspan, s0, options);

%% Offline observation 
% Constants 
Pobs = blkdiag( 1e-1 * eye(3), 1e-9 * eye(3), 0.1 );             % Covariance of the observation 

% Assume a linear identity observation model (H = I)
PD = 0.9;                                                        % Probability of observation
index = logical(randsrc(size(s,1), 1, [0, 1; 1-PD, PD]));        % Observation epochs

% Observed state vector
sobs = [tp(index) s(index,:)];

meas = sobs;
for i = 1:size(sobs,1)
    meas(i,2:end) = observation_model(sobs(i,2:end));   % Go through the measurement model
    meas(i,2:end) = mvnrnd(sobs(i,2:end), Pobs, 1);     % Add some noises
end

%% Offline estimation (EKF)
fprintf('Initialization... \n');
fprintf('------------------------------\n');
fprintf('Running the filter: \n');

% Constants 
Q = 1e-5 * eye(7);      % Dynamical noise 
R = 1e-4 * eye(7);      % Measurement noise 

s_dim = 7;              % Dimension of the state vector
y_dim = 7;              % Dimension of the observation vector

% Estimator setup
EKF_estimator = Filters.EKF();
EKF_estimator = EKF_estimator.AssignStateProcess(s_dim, @(Q, s, sigma, time_step)state_model(J2, Keci, Q, s, sigma, time_step));
EKF_estimator = EKF_estimator.AssignObservationProcess(y_dim, @(s)observation_model(s));
EKF_estimator = EKF_estimator.AdditiveCovariances(Q, R);

% Initial conditions 
se = [0; 0; 1.5; 0; 0; 0.01; 0];        % Estimator initial conditions 
Pe = 1e-2 * eye(7);                     % Initial covariance matrix

% Main loop
i = 1;
epoch = tp(1);                           % Initial epoch

while (i <= size(meas,1))
    % Propagation 
    tic 
    dt = meas(i,1) - epoch;
    if (dt > 0)
        EKF_estimator.State = se;
        EKF_estimator.Sigma = Pe;
        [se, Pe] = EKF_estimator.PropagationStep(dt);
    end

    % Correction
    [se, Pe] = EKF_estimator.CorrectionStep(se, Pe, meas(i,2:end).');

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
for i = 1:length(t)
    Sigma = reshape(P(:,i), s_dim, s_dim);
    [V, lambda] = eig(Sigma);
    
    sigma(:,i) = sqrt( diag(lambda) );
    e(:,i) = V * e(:,i);
end

e = e.';
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

subplot(3,1,2)
hold on
plot(tp, s(:,4:6))
scatter(meas(:,1), meas(:,5), 'Marker', 'x')
scatter(meas(:,1), meas(:,6), 'Marker', 'x')
scatter(meas(:,1), meas(:,7), 'Marker', 'x')
hold off
grid on
ylabel('$\|\mathbf{e}\|$ [-]')
legend('$e_x$', '$e_y$', '$e_z$')

subplot(3,1,3)
hold on
plot(tp, s(:,7), 'b')
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
plot(tp, sqrt( dot(s(:,4:6), s(:,4:6), 2) ), 'b')
ylabel('$\|\mathbf{e}\|$ [-]')
grid on
xlabel('$t$')

% Error plots 
figure 
subplot(3,1,1)
hold on
plot(t, e(:,1), 'b')
plot(t, e(:,1) + 3 * sigma(:,1))
plot(t, e(:,1) - 3 * sigma(:,1))
grid on
ylabel('$\Delta h_x$ [-]')
subplot(3,1,2)
hold on
plot(t, e(:,2), 'b')
plot(t, e(:,2) + 3 * sigma(:,2))
plot(t, e(:,2) - 3 * sigma(:,2))
grid on
ylabel('$\Delta h_y$ [-]')
subplot(3,1,3)
hold on
plot(t, e(:,3), 'b')
plot(t, e(:,3) + 3 * sigma(:,1))
plot(t, e(:,3) - 3 * sigma(:,1))
grid on
ylabel('$\Delta h_z$ [-]')
xlabel('$t$')

figure 
subplot(3,1,1)
plot(t, e(:,4), 'b')
plot(t, e(:,4) + 3 * sigma(:,4))
plot(t, e(:,4) - 3 * sigma(:,4))
grid on
ylabel('$\Delta e_x$ [-]')
subplot(3,1,2)
plot(t, e(:,5), 'b')
plot(t, e(:,5) + 3 * sigma(:,5))
plot(t, e(:,5) - 3 * sigma(:,5))
grid on
ylabel('$\Delta e_y$ [-]')
subplot(3,1,3)
plot(t, e(:,6), 'b')
plot(t, e(:,6) + 3 * sigma(:,6))
plot(t, e(:,6) - 3 * sigma(:,6))
grid on
ylabel('$\Delta e_z$ [-]')
xlabel('$t$')

figure 
hold on
plot(t, e(:,7), 'b')
plot(t, e(:,7) + 3 * sigma(:,7))
plot(t, e(:,7) - 3 * sigma(:,7))
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
function [s, P] = state_model(J2, Keci, Q, s0, P0, dt)
    % Integration of the state vector 
    options = odeset('AbsTol', 1E-22, 'RelTol', 2.25E-14);      % Integration tolerances

    [~, s] = ode45(@(t,s)Astrodynamics.milankovitch_dynamics(J2, Keci, t, s), [0 dt], s0, options);

    % Integration of the covariance matrix 
    A = Astrodynamics.milankovitch_jacobian(J2, Keci, s0);
    Phi = expm(A * dt);
    P = Phi * P0 * Phi.' + Q;
    s = s(end,:).';
end

% Observation model 
function [y, H] = observation_model(s)
    % Linear identity model
    y = s;                                   % Measurements
    H = eye(size(s,1), size(s,1));           % Observation matrix
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
