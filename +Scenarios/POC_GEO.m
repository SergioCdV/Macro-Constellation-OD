%% Combinatorial constellation macro-determination 
% Date: 02/02/2023
% Author: Sergio Cuevas del Valle

%% Constellation orbit determination. Scenario I %%
% This script provides a 1 plane constellation of 5 spacecraft, generating
% measurements from Madrid and Paris.

close all 
clear 

%% General user defined input
% Constants 
r0 = 6900e3;                % Characteristic distance of the Earth orbit    
mu = 3.86e14;               % Gravitional parameter of the Earth
Re = 6378e3;                % Reference Earth radius
J2 = 1e-3;                  % Earth's J2 parameter
Nmax = 1;                   % Number of targets

% Constellation lifetime
InitialEpoch = juliandate(datetime('now'));         % Initial epoch in JD
T = 2;                                              % Number of days 
EndEpoch = juliandate(datetime('now')+days(T));     % End epoch
Step = 600;                                         % Integration step in seconds
tspan = 0:Step:T * 86400;                           % Relative lifetime in seconds

% Target birth 
PS = 0.9999;                % Probability of surviving
PB = 0.00;                  % Birth rate
PD = 0.98;                  % Probability of detecting a target
Pc = 0.0;                   % Probability of false measurements
Vc = 10;                    % Number of false measurements per orbit 

%% Target births and deaths 
% Preallocation 
N = Nmax.*ones(1,length(tspan));    % Number of total spacecraft in time
S = cell(Nmax,1);                   % Time span of each target

% Original births
for i = 1:Nmax
    % Deaths of the original plane 
    deaths = logical(randsrc(length(tspan),1,[0, 1; PS, 1-PS]));
    td = find(deaths, 1, 'first'); 
      
    if (isempty(td))
        td = length(tspan);
    else
        N(1,td:end) = N(1,td:end)-1;
    end
    S{i} = [0 tspan(td)];  
end

GoOn = true;
i = 1;
while (GoOn)
    % Birth times  
    births = logical(randsrc(1,length(tspan)-2,[0, 1; 1-PB, PB])) & (N(1:end-2) < Nmax);
    tb = find(births, 1, 'first');

    % Death times
    if (~isempty(tb))
        tspan_d = logical(randsrc(1, length(tspan), [0, 1; PS, 1-PS])) & (tspan > tspan(tb+2));
        td = find(tspan_d, 1, 'first'); 
        if (isempty(td))
            td = length(tspan);
        end

        % Initial conditions
        S{i+Nmax} = [tspan(tb) tspan(td)];          
        i = i+1;
        N(1,tb:td) = N(1,tb:td)+1; 
    end

    % Convergence 
    if (all(N >= Nmax))
        GoOn = false;
    end
end

%% Constellation definition (only one plane) 
% Constellation definition 
Constellation_1 = Constellation('User defined', Nmax, 1, Nmax);
Constellation_1 = Constellation_1.ChangeTimeStep(Step);

% Orbit definition
ElementType = 'COE'; 
ElementSet = [r0 1e-3 0 deg2rad(0) deg2rad(0)]; 

for i = 1:size(S,1)
    % Generate a random anomaly 
    ElementSet = [ElementSet(1:5) 2*pi*rand()];

    % Add the orbit to the constellation
    AuxOrbit = Orbit(mu, ElementType, ElementSet, InitialEpoch + S{i}(1));
    AuxOrbit = AuxOrbit.SetFinalEpoch(EndEpoch + S{i}(2));
    AuxOrbit = AuxOrbit.Normalize(true, r0); 
    AuxOrbit = AuxOrbit.ChangeStateFormat('COE');
    AuxOrbit = AuxOrbit.DefineJ2Problem(J2, Re);
    AuxOrbit = AuxOrbit.AddPropagator('Osculating J2', Step);

    Constellation_1 = Constellation_1.AddOrbit(AuxOrbit);
end

% Set graphics
AuxOrbit.set_graphics();

% Compute the constellation parameters: number of planes, spacecraft and spacecraft per plane
Constellation_1.N = Constellation_1.NumberOfSpacecraft();
[Constellation_1.Np, Constellation_1.n] = Constellation_1.NumberOfPlanes();

% Propagation 
Constellation_1 = Constellation_1.Propagate(EndEpoch);

%% Observation process 
% Define the orbit to observer 
Orbit_2 = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('COE').Normalize(false, r0);

% Define an inertial observer
InObs = GibbsObserver().probability_detection(0.98).AddCovariance(diag([100 10])).AddInitialState(InitialEpoch, [1 0 0]);

% Prepare the measurements
[timestamp, meas, StateEvolution] = InObs.Observe(Orbit_2, EndEpoch);
Measurements = {[timestamp, meas], StateEvolution, @(meas, y)InObs.LikelihoodFunction(InObs.Sigma, meas, y), @(Orbit)InObs.ObservationProcess(timestamp, Orbit, StateEvolution)};

% Define a radar and telescope topocentric observer located at Madrid
RadarObs = TopocentricObserver('RADAR').probability_detection(0.98).AddFOV(deg2rad(120)).AddCovariance(diag([100 10])).AddInitialState(InitialEpoch, [deg2rad(40) deg2rad(-3)]);
TelescopeObs = TopocentricObserver('RADEC').probability_detection(0.98).AddFOV(deg2rad(120)).AddCovariance( deg2rad(diag([5 5])) ).AddInitialState(InitialEpoch, [deg2rad(40) deg2rad(-3)]);

% Prepare the measurements
[timestamp, meas_radar, StateEvolution] = RadarObs.Observe(Orbit_2, EndEpoch);
meas_radar = meas_radar + [normrnd(0,100,size(meas_radar,1),1) normrnd(0,10,size(meas_radar,1),1)];

% Clutter generation 
index = logical(randsrc(size(meas_radar,1), 1, [0, 1; 1-Pc, Pc]));
clutter = [normrnd(mean(meas_radar(:,1)), deg2rad(100), size(meas_radar, 1), 1) normrnd(mean(meas_radar(:,2)), deg2rad(10), size(meas_radar, 1), 1)]; 
clutter = [meas_radar(index,1) clutter(index,:)];
timestamp = [timestamp; timestamp(index,:)];
StateEvolution = [StateEvolution; StateEvolution(index,:)];

if (size(clutter,1) > Vc)
    index = randi([0 size(clutter,1)], Vc, 1); 
    clutter = clutter(index,:);
end

% Complete measurement set 
meas_radar = [meas_radar; clutter];
[~, index] = sort(meas_radar(:,1)); 
meas_radar = meas_radar(index,:);
StateEvolution = StateEvolution(index,:);
timestamp = timestamp(index);

Measurements_radar = {[timestamp, meas_radar], StateEvolution, @(meas, y)RadarObs.LikelihoodFunction(RadarObs.Sigma, meas, y), @(Orbit)RadarObs.ObservationProcess(timestamp, Orbit, StateEvolution)};

[timestamp, meas_radec, StateEvolution] = TelescopeObs.Observe(Orbit_2, EndEpoch);
meas_radec(:,1:2) = meas_radec(:,1:2) + [normrnd(0,deg2rad(5),size(meas_radec,1),1) normrnd(0,deg2rad(5),size(meas_radec,1),1)];

% Clutter generation 
index = logical(randsrc(size(meas_radec,1), 1, [0, 1; 1-Pc, Pc]));
clutter = normrnd(0, deg2rad(1), size(meas_radec, 1), 2);
clutter = [meas_radec(index,1) clutter(index,:)];
timestamp = [timestamp; timestamp(index,:)];
StateEvolution = [StateEvolution; StateEvolution(index,:)];

if (size(clutter,1) > Vc)
    index = randi([0 size(clutter,1)], Vc, 1); 
    clutter = clutter(index,:);
end

% Complete measurement set 
meas_radec = [meas_radec; clutter];
[~, index] = sort(meas_radec(:,1)); 
meas_radec = meas_radec(index,:);
StateEvolution = StateEvolution(index,:);
timestamp = timestamp(index);

Measurements_telescope = {[timestamp, meas_radec], StateEvolution, @(meas, y)TelescopeObs.LikelihoodFunction(RadarObs.Sigma, meas, y), @(Orbit)TelescopeObs.ObservationProcess(timestamp, Orbit, StateEvolution)};

%% Estimation 
% Estimator configuration

% Estimation

%% Results


%% Auxiliary functions
% Kinematic proposal 
function [theta_plus, sigma_plus] = kinematic_proposal(mu, r, sigma_r, time_step, theta, sigma)
    % State update
    theta_plus = theta + time_step*sqrt(mu/r^3);

    % Covariance update
    sigma_plus = sigma + (9/4)*time_step^2*r^(-5)*sigma_r;
end

% Observation proess 
function [y, H] = radar(theta)
    % Expected observation
    y = [cos(theta); sin(theta)];

    % Observation matrix 
    H = [-sin(theta); cos(theta)];
end

% Compute the likelihood function 
function [q] = likelihood_function(z, m, P)
    q = exp(-0.5*(z-m).'*P^(-1)*(z-m))/sqrt(det(P)*(2*pi)^size(P,1));
end