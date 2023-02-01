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
EndEpoch = juliandate(datetime('tomorrow'));

% Orbit definition
ElementType = 'COE'; 
ElementSet = [r0 1e-3 0 deg2rad(0) deg2rad(0) deg2rad(0)]; 

Orbit = Orbit(mu, ElementType, ElementSet, InitialEpoch);
Orbit = Orbit.Normalize(true, r0);
Orbit = Orbit.Normalize(true, 2*r0);
Orbit = Orbit.SetFinalEpoch(EndEpoch); 

% Transformation of elements 
Orbit = Orbit.ChangeStateFormat('MOE');
Orbit = Orbit.ChangeStateFormat('COE');
Orbit = Orbit.ChangeStateFormat('Cartesian');
Orbit = Orbit.ChangeStateFormat('COE');
Orbit = Orbit.ChangeStateFormat('Cartesian');
Orbit = Orbit.ChangeStateFormat('KS');
Orbit = Orbit.ChangeStateFormat('MOE');
Orbit = Orbit.ChangeStateFormat('KS');
Orbit = Orbit.ChangeStateFormat('Cartesian');
Orbit = Orbit.ChangeStateFormat('MOE');
Orbit = Orbit.ChangeStateFormat('KS');
Orbit = Orbit.ChangeStateFormat('MOE');
Orbit = Orbit.ChangeStateFormat('COE');
Orbit = Orbit.ChangeStateFormat('KS');
Orbit = Orbit.ChangeStateFormat('Cartesian');

% Orbit propagation 
Orbit = Orbit.AddPropagator('Keplerian', 1E-1);
Orbit = Orbit.SetCurrentEpoch(juliandate(datetime('now'))+5000);
Orbit = Orbit.Propagate();

% Set graphics 
Orbit.set_graphics();

% Set trajectory 
Orbit.PlotTrajectory(figure(1), Orbit.InitialEpoch, Orbit.CurrentEpoch)