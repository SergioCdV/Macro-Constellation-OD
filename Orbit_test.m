%% Combinatorial constellation macro-determination 
% Date: 31/01/2022
% Author: Sergio Cuevas del Valle

%% Orbit test %% 
% Script to test the functionalities of the orbit class

close all; 
clear;

%% Orbit definition 
% Constants of the environment 
r0 = 149597870700;                      % 1 AU [m]
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2] 

r0 = 6900e3; 
mu = 3.86e14;

% Epochs
InitialEpoch = juliandate(datetime('now'));
EndEpoch = juliandate(datetime('tomorrow')+60000);

% Orbit definition
ElementType = 'COE'; 
ElementSet = [r0 1e-3 0 deg2rad(90) deg2rad(0) deg2rad(0)]; 

Orbit = Orbit(mu, ElementType, ElementSet, InitialEpoch);
% Orbit = Orbit.Normalize(true, r0);
% Orbit = Orbit.Normalize(true, 2*r0);
Orbit = Orbit.SetFinalEpoch(EndEpoch); 

% Transformation of elements 
Orbit = Orbit.ChangeStateFormat('MOE');
Orbit = Orbit.ChangeStateFormat('COE');
Orbit = Orbit.ChangeStateFormat('Cartesian');
Orbit = Orbit.ChangeStateFormat('COE');
Orbit = Orbit.ChangeStateFormat('Cartesian');
Orbit = Orbit.ChangeStateFormat('KS');
Orbit = Orbit.ChangeStateFormat('MOE');
Orbit = Orbit.ChangeStateFormat('COE');
Orbit = Orbit.ChangeStateFormat('KS');
Orbit = Orbit.ChangeStateFormat('COE');
Orbit = Orbit.ChangeStateFormat('Cartesian');
Orbit = Orbit.ChangeStateFormat('MOE');
Orbit = Orbit.ChangeStateFormat('KS');
Orbit = Orbit.ChangeStateFormat('MOE');
Orbit = Orbit.ChangeStateFormat('COE');
Orbit = Orbit.ChangeStateFormat('KS');
Orbit = Orbit.ChangeStateFormat('Cartesian');

% Orbit propagation 
Orbit = Orbit.AddPropagator('Keplerian', 0.5);

% Set trajectory 
Orbit.set_graphics();

%% Constellation definition
% Constellation constructor
Constellation_test = Constellation('User defined', 4,4,4);

% Add/remove orbits
Constellation_test = Constellation_test.AddOrbit(Orbit);
Constellation_test = Constellation_test.AddOrbit(Orbit);
Constellation_test = Constellation_test.AddOrbit(Orbit);
Constellation_test = Constellation_test.AddOrbit(Orbit);

Orbit_2 = Orbit; 
Orbit_2.ElementSet(3) = pi; 
Orbit_2.SetCurrentEpoch(Orbit_2.InitialEpoch); 
Constellation_test = Constellation_test.AddOrbit(Orbit_2);

% Compute number of planes, spacecraft and spacecraft per plane
Constellation_test.N = Constellation_test.NumberOfSpacecraft();
[Constellation_test.Np, Constellation_test.n] = Constellation_test.NumberOfPlanes();

% Propagation 
Constellation_test = Constellation_test.Propagate(juliandate(datetime('now'))+6000);
Constellation_test.OrbitSet{1,2}.PlotTrajectory(figure(1), Constellation_test.OrbitSet{1,2}.InitialEpoch, Constellation_test.OrbitSet{1,2}.FinalEpoch);
