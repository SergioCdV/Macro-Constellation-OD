%% Combinatorial constellation macro-determination 
% Date: 31/01/2022
% Author: Sergio Cuevas del Valle

%% Orbit test %% 
% Script to test the functionalities of the fundamental problem classes

close all; 
clear;

%% Orbit definition 
% Constants of the environment 
r0 = 6900e3;                  % 1 AU [m]
mu = 3.86e14;                 % Gavitational parameter of the Sun [m^3 s^−2] 

% Epochs
InitialEpoch = juliandate(datetime('now'));
EndEpoch = juliandate(datetime('tomorrow') + days(7));

% Orbit definition
ElementType = 'COE'; 
ElementSet = [r0 1e-3 deg2rad(10) deg2rad(5) deg2rad(10) deg2rad(0)]; 

% ElementSet = Astrodynamics.sso_elements(6, ElementSet);
% ElementSet(2) = 1e-2;
% ElementSet = [ElementSet ElementSet(1) * sqrt(1-ElementSet(2)^2)];

Orbit_1 = Orbit(mu, ElementType, ElementSet, InitialEpoch).Normalize(true, r0);
Orbit_1 = Orbit_1.SetFinalEpoch(EndEpoch); 

% Transformation of elements 
% Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
% Orbit_1 = Orbit_1.ChangeStateFormat('COE');
% Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
% Orbit_1 = Orbit_1.ChangeStateFormat('COE');
% Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
% Orbit_1 = Orbit_1.ChangeStateFormat('KS');
% Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
% Orbit_1 = Orbit_1.ChangeStateFormat('COE');
% Orbit_1 = Orbit_1.ChangeStateFormat('KS');
% Orbit_1 = Orbit_1.ChangeStateFormat('COE');
% Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
% Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
% Orbit_1 = Orbit_1.ChangeStateFormat('KS');
% Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
% Orbit_1 = Orbit_1.ChangeStateFormat('COE');
% Orbit_1 = Orbit_1.ChangeStateFormat('KS');
% Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
% Orbit_1 = Orbit_1.ChangeStateFormat('POL');
% Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
% Orbit_1 = Orbit_1.ChangeStateFormat('POL');
% Orbit_1 = Orbit_1.ChangeStateFormat('KS');
% Orbit_1 = Orbit_1.ChangeStateFormat('POL');
% Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
% Orbit_1 = Orbit_1.ChangeStateFormat('POL');
% Orbit_1 = Orbit_1.ChangeStateFormat('COE');
% Orbit_1 = Orbit_1.ChangeStateFormat('POL');
% Orbit_1 = Orbit_1.ChangeStateFormat('COE');

% Orbit propagation 
Orbit_1 = Orbit_1.AddPropagator('Keplerian', 1);

% Set trajectory 
Orbit_1.set_graphics();

%% J2 problem testing
% Constants of the environment 
J2 = 1.08263e-3;            % Second zonal harmonic of the Earth
Re = 6378.14e3;             % Reference J2 radius [m] 

% Orbit definition
ElementType = 'COE'; 
Orbit_2 = Orbit(mu, ElementType, ElementSet, InitialEpoch);

Orbit_2 = Orbit_2.DefineJ2Problem(J2, Re);
Orbit_2 = Orbit_2.SetCurrentEpoch(EndEpoch); 
Orbit_2 = Orbit_2.Normalize(true, r0);
Orbit_3 = Orbit_2;
Orbit_4 = Orbit_3;

% Orbit propagation 
Orbit_2 = Orbit_2.AddPropagator('Mean J2', 60);
Orbit_3 = Orbit_3.AddPropagator('Osculating J2', 60);
Orbit_4 = Orbit_4.AddPropagator('High-precision', 1);


Orbit_2 = Orbit_2.Propagate();
Orbit_3 = Orbit_3.Propagate();
% Orbit_4 = Orbit_4.Propagate();

for j = 1:size(Orbit_2.StateEvolution,1)
    a = Orbit_2.StateEvolution(j,2);
    e = Orbit_2.StateEvolution(j,3);
    Omega = Orbit_2.StateEvolution(j,4);
    i = Orbit_2.StateEvolution(j,5);
    omega = Orbit_2.StateEvolution(j,6);
    M = Orbit_2.StateEvolution(j,7);
    D = [M omega Omega sqrt(a) sqrt(a*(1-e^2)) sqrt(a*(1-e^2))*cos(i)];
    Lara_orbit(:,j) = Astrodynamics.Lara_solution(J2, D);
end
%%
hold on;
% Orbit_2.PlotTrajectory(figure(1), Orbit_2.InitialEpoch, Orbit_2.PropagatedEpoch);
Orbit_3.PlotTrajectory(figure(1), Orbit_3.InitialEpoch, Orbit_3.PropagatedEpoch);
plot3(Lara_orbit(1,:), Lara_orbit(2,:), Lara_orbit(3,:))
% Orbit_4.PlotTrajectory(figure(1), Orbit_4.InitialEpoch, Orbit_4.PropagatedEpoch);
hold off
%%
figure
diff = Orbit_3.ChangeStateFormat('ECI').StateEvolution(:,2:end)-Orbit_2.ChangeStateFormat('ECI').StateEvolution(:,2:end);
diff2 = Orbit_3.ChangeStateFormat('ECI').StateEvolution(:,2:end)-Lara_orbit.';
hold on 
plot(Orbit_2.StateEvolution(:,1), log(sqrt(dot(diff,diff,2))))
plot(Orbit_2.StateEvolution(:,1), log(sqrt(dot(diff2,diff2,2))))
hold off 
legend('mean', 'Lara')

%% Results 
figure 
hold on
plot(Orbit_2.StateEvolution(:,1), Orbit_2.StateEvolution(:,4));
plot(Orbit_3.StateEvolution(:,1), Orbit_3.StateEvolution(:,4));
hold off
xlabel('$t$')
ylabel('$\Omega$')
legend('$\hat{\Omega}$', '$\Omega$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));
%%
figure 
hold on
plot(Orbit_2.StateEvolution(:,1), Orbit_2.StateEvolution(:,6), 'k', 'Linewidth', 1);
plot(Orbit_3.StateEvolution(1:size(Orbit_3.StateEvolution,1),1), Orbit_3.StateEvolution(1:size(Orbit_3.StateEvolution,1),6), 'Linewidth', 0.1);
hold off
xlabel('$t$')
ylabel('$\omega$')
legend('$\hat{\omega}$', '$\omega$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));
%%
figure 
hold on
plot(Orbit_2.StateEvolution(:,1), log(Orbit_2.StateEvolution(:,4)-Orbit_3.StateEvolution(:,4)));
hold off
xlabel('$t$')
ylabel('log $e_{\Omega}$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));

figure 
hold on
plot(Orbit_2.StateEvolution(:,1), log(Orbit_2.StateEvolution(:,6)-Orbit_3.StateEvolution(:,6)));
hold off
xlabel('$t$')
ylabel('log $e_{\omega}$')
grid on;
yticklabels(strrep(yticklabels, '-', '$-$'));

%% Constellation definition
% Constellation constructor
Constellation_test = Constellation('User defined', 4,4,4);

% Add/remove orbits
Constellation_test = Constellation_test.AddOrbit(Orbit_1);
Constellation_test = Constellation_test.AddOrbit(Orbit_1);
Constellation_test = Constellation_test.AddOrbit(Orbit_1);
Constellation_test = Constellation_test.AddOrbit(Orbit_1);

Orbit_2 = Orbit_1; 
Orbit_2.ElementSet(3) = pi; 
Orbit_2.SetCurrentEpoch(Orbit_2.InitialEpoch); 
Constellation_test = Constellation_test.AddOrbit(Orbit_2);

% Compute number of planes, spacecraft and spacecraft per plane
Constellation_test.N = Constellation_test.NumberOfSpacecraft();
[Constellation_test.Np, Constellation_test.n] = Constellation_test.NumberOfPlanes();

% Propagation 
Constellation_test = Constellation_test.Propagate(juliandate(datetime('now'))+6000);
Constellation_test.OrbitSet{1,2}.PlotTrajectory(figure(1), Constellation_test.OrbitSet{1,2}.InitialEpoch, Constellation_test.OrbitSet{1,2}.FinalEpoch);

%% Observations 
% Define an inertial observer
InitialState = [0 1 0];
InitialEpoch = juliandate(datetime('now'));
Sigma = diag([1e3 1e3 1e3]);
PD = 0.98;
InObs = Sensors.GibbsSensor(InitialEpoch, InitialState, Sigma, PD);

% Prepare the measurements
FinalObserveEpoch = juliandate(datetime('now'))+1800/(24 * 3600);
[timestamp, meas, StateEvolution] = InObs.Observe(Orbit_2, FinalObserveEpoch);

Measurements = {[timestamp, meas], StateEvolution, @(meas, y)InObs.LikelihoodFunction(InObs.Sigma, meas, y), @(Orbit)InObs.ObservationProcess(timestamp, Orbit, StateEvolution)};

% Define a radar topocentric observer located at Madrid
InitialState = [deg2rad(40) deg2rad(-3)];
InitialEpoch = juliandate(datetime('now'));
Sigma = eye(2,2);
PD = 0.98;
FOV = deg2rad(120);
RadarObs = Sensors.RadarSensor(InitialEpoch, InitialState, Sigma, PD, FOV);

% Prepare the measurements
FinalObserveEpoch = juliandate(datetime('now'))+1800/(24 * 3600);
[timestamp, meas_radar, StateEvolution] = RadarObs.Observe(Orbit_2, FinalObserveEpoch);

Measurements_radar = {[timestamp, meas_radar], StateEvolution, @(meas, y)RadarObs.LikelihoodFunction(RadarObs.Sigma, meas, y), @(Orbit)RadarObs.ObservationProcess(timestamp, Orbit, StateEvolution)};

% Define a telescope topocentric observer located at Madrid
InitialState = [deg2rad(40) deg2rad(-3)];
InitialEpoch = juliandate(datetime('now'));
Sigma = diag([deg2rad(1) deg2rad(1) deg2rad(0.00001) deg2rad(0.00001)]);
PD = 0.98;
FOV = deg2rad(120);
TelescopeObs = Sensors.TopocentricSensor(InitialEpoch, InitialState, Sigma, PD);

% Prepare the measurements
FinalObserveEpoch = juliandate(datetime('now'))+1800/(24 * 3600);
[timestamp, meas_radec, StateEvolution] = TelescopeObs.Observe(Orbit_2, FinalObserveEpoch);

Measurements_telescope = {[timestamp, meas_radec], StateEvolution, @(meas, y)TelescopeObs.LikelihoodFunction(RadarObs.Sigma, meas, y), @(Orbit)TelescopeObs.ObservationProcess(timestamp, Orbit, StateEvolution)};