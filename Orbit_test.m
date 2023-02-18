%% Combinatorial constellation macro-determination 
% Date: 31/01/2022
% Author: Sergio Cuevas del Valle

%% Orbit test %% 
% Script to test the functionalities of the orbit class

close all; 
clear;

setup_path();

%% Orbit definition 
% Constants of the environment 
r0 = 149597870700;                      % 1 AU [m]
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2] 

r0 = 6900e3; 
mu = 3.86e14;

% Epochs
InitialEpoch = juliandate(datetime('now'));
EndEpoch = juliandate(datetime('tomorrow'));

% Orbit definition
ElementType = 'COE'; 
ElementSet = [r0 1e-3 0 deg2rad(90) deg2rad(0) deg2rad(0)]; 

Orbit_1 = Orbit(mu, ElementType, ElementSet, InitialEpoch).Normalize(true, 2*r0);
Orbit_1 = Orbit_1.SetFinalEpoch(EndEpoch); 

% Transformation of elements 
Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
Orbit_1 = Orbit_1.ChangeStateFormat('COE');
Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
Orbit_1 = Orbit_1.ChangeStateFormat('COE');
Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
Orbit_1 = Orbit_1.ChangeStateFormat('KS');
Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
Orbit_1 = Orbit_1.ChangeStateFormat('COE');
Orbit_1 = Orbit_1.ChangeStateFormat('KS');
Orbit_1 = Orbit_1.ChangeStateFormat('COE');
Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
Orbit_1 = Orbit_1.ChangeStateFormat('KS');
Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
Orbit_1 = Orbit_1.ChangeStateFormat('COE');
Orbit_1 = Orbit_1.ChangeStateFormat('KS');
Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
Orbit_1 = Orbit_1.ChangeStateFormat('POL');
Orbit_1 = Orbit_1.ChangeStateFormat('ECI');
Orbit_1 = Orbit_1.ChangeStateFormat('POL');
Orbit_1 = Orbit_1.ChangeStateFormat('KS');
Orbit_1 = Orbit_1.ChangeStateFormat('POL');
Orbit_1 = Orbit_1.ChangeStateFormat('MOE');
Orbit_1 = Orbit_1.ChangeStateFormat('POL');
Orbit_1 = Orbit_1.ChangeStateFormat('COE');
Orbit_1 = Orbit_1.ChangeStateFormat('POL');
Orbit_1 = Orbit_1.ChangeStateFormat('COE');

% Orbit propagation 
Orbit_1 = Orbit_1.AddPropagator('Keplerian', 0.5);

% Set trajectory 
Orbit_1.set_graphics();

%% J2 problem testing
% Constants of the environment 
J2 = 1.08263e-3;            % Second zonal harmonic of the Earth
Re = 6378.14e3;             % Reference J2 radius [m] 

% Orbit definition
ElementType = 'COE'; 
ElementSet = [r0 1e-3 0 deg2rad(90) deg2rad(0) deg2rad(0)]; 

Orbit_2 = Orbit(mu, ElementType, ElementSet, InitialEpoch);

Orbit_2 = Orbit_2.DefineJ2Problem(J2, Re);
Orbit_2 = Orbit_2.Normalize(true, r0);
Orbit_2 = Orbit_2.SetFinalEpoch(EndEpoch); 
Orbit_3 = Orbit_2;
Orbit_4 = Orbit_3;

% Orbit propagation 
Orbit_2 = Orbit_2.AddPropagator('Mean J2', 0.5);
Orbit_3 = Orbit_3.AddPropagator('Osculating J2', 0.5);
Orbit_3 = Orbit_3.AddPropagator('High-precision', 0.5);

Orbit_2 = Orbit_2.SetCurrentEpoch(EndEpoch);
Orbit_2 = Orbit_2.Propagate();
Orbit_2.set_graphics();
Orbit_2.PlotTrajectory(figure(1), Orbit_2.InitialEpoch, Orbit_2.PropagatedEpoch);

Orbit_3 = Orbit_3.SetCurrentEpoch(EndEpoch);
Orbit_3 = Orbit_3.Propagate();

Orbit_4 = Orbit_4.SetCurrentEpoch(EndEpoch);
Orbit_4 = Orbit_4.Propagate();

Orbit_3.set_graphics();
hold on;
Orbit_3.PlotTrajectory(figure(1), Orbit_3.InitialEpoch, Orbit_3.PropagatedEpoch);
hold on;
Orbit_4.PlotTrajectory(figure(1), Orbit_4.InitialEpoch, Orbit_4.PropagatedEpoch);

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
InObs = GibbsObserver().probability_detection(0.98).AddCovariance(eye(2,2)).AddInitialState(0, [1 0 0]);

% Prepare the measurements
[timestamp, meas, StateEvolution] = InObs.Observe(Orbit_2, juliandate(datetime('now')));

Measurements = {[timestamp, meas], StateEvolution, @(meas, y)InObs.LikelihoodFunction(InObs.Sigma, meas, y), @(Orbit)InObs.ObservationProcess(timestamp, Orbit, StateEvolution)};

% Define a radar topocentric observer located at Madrid
RadarObs = TopocentricObserver('RADAR').probability_detection(0.98).AddFOV(deg2rad(120)).AddCovariance(eye(2,2)).AddInitialState(juliandate(datetime('now')), [deg2rad(40) deg2rad(-3)]);

% Prepare the measurements
[timestamp, meas_radar, StateEvolution] = RadarObs.Observe(Orbit_2, juliandate(datetime('now'))+1800/(24 * 3600));

Measurements_radar = {[timestamp, meas_radar], StateEvolution, @(meas, y)RadarObs.LikelihoodFunction(RadarObs.Sigma, meas, y), @(Orbit)RadarObs.ObservationProcess(timestamp, Orbit, StateEvolution)};

% Define a telescope topocentric observer located at Madrid
TelescopeObs = TopocentricObserver('RADEC').probability_detection(0.98).AddFOV(deg2rad(120)).AddCovariance(eye(2,2)).AddInitialState(juliandate(datetime('now')), [deg2rad(40) deg2rad(-3)]);

% Prepare the measurements
[timestamp, meas_radec, StateEvolution] = TelescopeObs.Observe(Orbit_2, juliandate(datetime('now'))+1800/(24 * 3600));

Measurements_telescope = {[timestamp, meas_radec], StateEvolution, @(meas, y)TelescopeObs.LikelihoodFunction(RadarObs.Sigma, meas, y), @(Orbit)TelescopeObs.ObservationProcess(timestamp, Orbit, StateEvolution)};