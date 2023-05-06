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
Nmax = 12;                  % Number of targets

% Constellation lifetime
InitialEpoch = juliandate(datetime('now'));         % Initial epoch in JD
T = 2;                                              % Number of days 
EndEpoch = juliandate(datetime('now')+days(T));     % End epoch
Step = 600;                                         % Integration step in seconds
tspan = 0:Step:T * 86400;                           % Relative lifetime in seconds

% Target birth 
PS = 1;                     % Probability of surviving
PB = 0.00;                  % Birth rate 

%% Target births and deaths 
% Preallocation 
N = Nmax.*ones(1,length(tspan));    % Number of total spacecraft in time
S = cell(Nmax,1);                   % Time span of each target

% Original births
for i = 1:Nmax
    % Deaths of the original set 
    deaths = logical(randsrc(length(tspan), 1, [0, 1; PS, 1-PS]));
    td = find(deaths, 1, 'first'); 
      
    if (isempty(td))
        td = length(tspan);
    else
        N(1,td:end) = max(0, N(1,td:end)-1);
    end 
    S{i} = [0 max(tspan(td), Step)];
end

GoOn = true;
i = 1;
while (GoOn)
    % Birth times  
    births = logical(randsrc(1, length(tspan)-1, [0, 1; 1-PB, PB])) & (N(1:end-1) < Nmax);
    tb = find(births, 1, 'first');

    % Death times
    if (~isempty(tb))
        tspan_d = logical(randsrc(1, length(tspan(tb:end)), [0, 1; PS, 1-PS]));
        td = find(tspan_d, 1, 'first'); 
        if (isempty(td))
            td = length(tspan);
        else
            td = max(td,2);
            td = min(tb+td, length(tspan));
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
Planes(1,:) = [r0 1e-3 deg2rad(0) deg2rad(30) deg2rad(270)]; 
Planes(2,:) = [r0 1e-3 deg2rad(120) deg2rad(90) deg2rad(270)]; 
Planes(3,:) = [r0 1e-3 deg2rad(240) deg2rad(150) deg2rad(270)]; 

ElementSet = [Scenarios.sso_elements(6, repmat(Planes(1, [1 2 5]),  4, 1)); ...
              Scenarios.sso_elements(10, repmat(Planes(1, [1 2 5]), 4, 1));
              Scenarios.sso_elements(12, repmat(Planes(1, [1 2 5]), 4, 1))];

for i = 1:size(ElementSet,1)
    % Add the orbit to the constellation
    AuxOrbit = Orbit(mu, ElementType, ElementSet(i,:), InitialEpoch + S{i}(1)/86400);
    AuxOrbit = AuxOrbit.SetFinalEpoch(InitialEpoch + S{i}(2)/86400);
    AuxOrbit = AuxOrbit.Normalize(true, r0); 
    AuxOrbit = AuxOrbit.ChangeStateFormat('COE');
    AuxOrbit = AuxOrbit.DefineJ2Problem(J2, Re);
    AuxOrbit = AuxOrbit.AddPropagator('Osculating J2', Step);

    Constellation_1 = Constellation_1.AddOrbit(AuxOrbit);
end

%% Sensor network 
% Clutter
Pc = 0.0;                   % Probability of false measurements
Vc = 10;                    % Number of false measurements per sensor (surveillance region)

% Define a radar topocentric observer located at Madrid
InitialState = [deg2rad(40) deg2rad(-3)];
InitialEpoch = juliandate(datetime('now'));
Sigma = eye(2,2);
PD = 0.98;
FOV = deg2rad(120);
RadarObs = Sensors.RadarSensor(InitialEpoch, InitialState, Sigma, PD, FOV);

% Define a telescope topocentric observer located at Madrid
InitialState = [deg2rad(40) deg2rad(-3)];
InitialEpoch = juliandate(datetime('now'));
Sigma = diag([deg2rad(1) deg2rad(1) deg2rad(0.00001) deg2rad(0.00001)]);
PD = 0.98;
FOV = deg2rad(120);
TelescopeObs_1 = Sensors.TopocentricSensor(InitialEpoch, InitialState, Sigma, PD);

% Define a telescope topocentric observer located at Madrid
InitialState = [deg2rad(40) deg2rad(-3)];
InitialEpoch = juliandate(datetime('now'));
Sigma = diag([deg2rad(1) deg2rad(1) deg2rad(0.00001) deg2rad(0.00001)]);
PD = 0.98;
FOV = deg2rad(120);
TelescopeObs_2 = Sensors.TopocentricSensor(InitialEpoch, InitialState, Sigma, PD);

% Define a telescope topocentric observer located at Madrid
InitialState = [deg2rad(40) deg2rad(-3)];
InitialEpoch = juliandate(datetime('now'));
Sigma = diag([deg2rad(1) deg2rad(1) deg2rad(0.00001) deg2rad(0.00001)]);
PD = 0.98;
FOV = deg2rad(120);
TelescopeObs_3 = Sensors.TopocentricSensor(InitialEpoch, InitialState, Sigma, PD);

%% Observation process 
% Prepare the measurements
FinalObserveEpoch = EndEpoch;

RadarTime = [];
meas_radar = [];
RadarState = [];

TelescopeTime_1 = [];
meas_radec_1 = [];
TelescopeState_1 = [];

TelescopeTime_2 = [];
meas_radec_2 = [];
TelescopeState_2 = [];

TelescopeTime_3 = [];
meas_radec_3 = [];
TelescopeState_3 = [];

for i = 1:size(Constellation_1.OrbitSet,1)
    [RadarTime_aux, meas_radar_aux, RadarState_aux] = RadarObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

    RadarTime = [RadarTime; RadarTime_aux];
    meas_radar = [meas_radar; meas_radar_aux];
    RadarState = [RadarState; RadarState_aux];

    [TelescopeTime_aux, meas_radec_aux, TelescopeState_aux] = TelescopeObs_1.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

    TelescopeTime_1 = [TelescopeTime_1; TelescopeTime_aux];
    meas_radec_1 = [meas_radec_1; meas_radec_aux];
    TelescopeState_1 = [TelescopeState_1; TelescopeState_aux];

    [TelescopeTime_aux, meas_radec_aux, TelescopeState_aux] = TelescopeObs_2.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

    TelescopeTime_2 = [TelescopeTime_2; TelescopeTime_aux];
    meas_radec_2 = [meas_radec_2; meas_radec_aux];
    TelescopeState_2 = [TelescopeState_2; TelescopeState_aux];

    [TelescopeTime_aux, meas_radec_aux, TelescopeState_aux] = TelescopeObs_3.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

    TelescopeTime_3 = [TelescopeTime_3; TelescopeTime_aux];
    meas_radec_3 = [meas_radec_3; meas_radec_aux];
    TelescopeState_3 = [TelescopeState_3; TelescopeState_aux];
end

ObservationSpan = [RadarTime; TelescopeTime_1; TelescopeTime_2; TelescopeTime_3];
[ObservationSpan, index] = sort(ObservationSpan);

Measurements = cell(length(ObservationSpan), 4);
Measurements(:,1) = num2cell(ObservationSpan);

for i = 1:length(index)
    if (index(i) <= size(RadarTime,1))
        Measurements(i,2) = { meas_radar(index(i),:) };
        Measurements(i,3) = { RadarState(index(i),:) };
        Measurements(i,4) = { @(y)RadarObs.LikelihoodFunction(RadarObs.Sigma, meas_radar(index(i),:), y) };

    elseif (index(i) <= size(RadarTime,1) + size(TelescopeTime_1,1))
        L = size(RadarTime,1);
        Measurements(i,2) = { meas_radec_1(index(i)-L,:) };
        Measurements(i,3) = { TelescopeState_1(index(i)-L,:) };
        Measurements(i,4) = { @(y)TelescopObs_1.LikelihoodFunction(TelescopObs_1.Sigma, meas_radec_1(index(i)-L,:), y) };

    elseif (index(i) <= size(RadarTime,1) + size(TelescopeTime_1,1) + size(TelescopeTime_2,1))
        L = size(RadarTime,1) + size(TelescopeTime_1,1);
        Measurements(i,2) = { meas_radec_2(index(i)-L,:) };
        Measurements(i,3) = { TelescopeState_2(index(i)-L,:) };
        Measurements(i,4) = { @(y)TelescopObs_2.LikelihoodFunction(TelescopObs_2.Sigma, meas_radec_2(index(i)-L,:), y) };

    else
        L = size(RadarTime,1) + size(TelescopeTime_1,1) + size(TelescopeTime_2,1);
        Measurements(i,2) = { meas_radec_3(index(i)-L,:) };
        Measurements(i,3) = { TelescopeState_2(index(i)-L,:) };
        Measurements(i,4) = { @(y)TelescopeObs_3.LikelihoodFunction(TelescopeObs_3.Sigma, meas_radec(index(i)-L,:), y) };
    end
end

%% Propagation 
% Propagation 
Constellation_1 = Constellation_1.Propagate(EndEpoch);

% Set graphics
AuxOrbit.set_graphics();

% Compute the constellation parameters: number of planes, spacecraft and spacecraft per plane
Constellation_1.N = Constellation_1.NumberOfSpacecraft();
[Constellation_1.Np, Constellation_1.n] = Constellation_1.NumberOfPlanes();

%% Estimation 
% Estimator configuration

% Estimation

%% Analysis

%% Results
figure
hold on
plot(tspan, Constellation_1.N)
%scatter(tspan, ); 
hold off
legend('$N$',' $\hat{N}$');
xlabel('Epoch $t$')
ylabel('N. SC')
grid on;

figure
hold on
plot(tspan, Constellation_1.Np)
%scatter(tspan, ); 
hold off
legend('$N_p$',' $\hat{N}_p$');
xlabel('Epoch $t$')
ylabel('N. planes')
grid on;

figure
hold on
plot(tspan, Constellation_1.n(:,1))
%scatter(tspan, ); 
legend('$n_{i,p}$',' $\hat{n}_{i,p}$', 'AutoUpdate', 'off');
if (size(Constellation_1.n,2) > 1)
    plot(tspan, Constellation_1.n(:,2:end))
    %scatter(tspan, )
end
hold off
xlabel('Epoch $t$')
ylabel('N. SC')
grid on;
%%
% Constellation plot 
figure(10)
hold on

Earth().plot();

Constellation_1.OrbitSet{1,2}.Normalize(false,r0).PlotTrajectory(figure(10), InitialEpoch, EndEpoch);
Constellation_1.OrbitSet{5,2}.Normalize(false,r0).PlotTrajectory(figure(10), InitialEpoch, EndEpoch);
Constellation_1.OrbitSet{9,2}.Normalize(false,r0).PlotTrajectory(figure(10), InitialEpoch, EndEpoch);

for i = 1:size(Constellation_1.OrbitSet,1)
    Raux = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('ECI').Normalize(false,r0).StateEvolution(end,1:4);
    if (Raux(1) >= tspan(end))
        Raux = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('ECI').Normalize(false,r0).StateEvolution(end,2:4);
        scatter3(Raux(1), Raux(2), Raux(3), 50, 'filled', 'r', 'LineWidth', 2);
    end
end
hold off
xlabel('$X$')
ylabel('$Y$')
zlabel('$Z$')
grid on;

% Anomaly plot 

% Distribution on the sphere (attitude planes)

% Distribution on the momenta space

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