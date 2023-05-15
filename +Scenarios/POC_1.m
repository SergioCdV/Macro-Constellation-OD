%% Combinatorial constellation macro-determination 
% Date: 02/02/2023
% Author: Sergio Cuevas del Valle

%% Constellation orbit determination. Scenario I %%
% This script provides a 1 plane constellation of 5 spacecraft, generating
% measurements from Madrid and Paris.

close all 
clear 
rng(1);

%% General user defined input
% Constants 
r0 = 6900e3;                % Characteristic distance of the Earth orbit    
mu = 3.986e14;              % Gravitional parameter of the Earth
Re = 6378e3;                % Reference Earth radius
Tc = sqrt(Re^2/mu);         % Characteristic time
J2 = 1e-3;                  % Earth's J2 parameter
Nmax = 4;                   % Number of targets

% Constellation lifetime
InitialEpoch = juliandate(datetime('now'));         % Initial epoch in JD
T = 1;                                              % Number of days 
EndEpoch = juliandate(datetime('now')+days(T));     % End epoch
Step = 600;                                          % Integration step in seconds
tspan = 0:Step:T * 86400;                           % Relative lifetime in seconds

% Target birth 
PS = 0.99;                  % Probability of surviving
PB = 0.01;                  % Birth rate 

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
ElementSet = [r0 1e-3 deg2rad(0) deg2rad(45) deg2rad(0)]; 

for i = 1:size(S,1)
    % Generate a random anomaly 
    ElementSet = [ElementSet(1:5) 2*pi*rand()];

    % Add the orbit to the constellation
    AuxOrbit = Orbit(mu, ElementType, ElementSet, InitialEpoch + S{i}(1)/86400);
    AuxOrbit = AuxOrbit.SetFinalEpoch(InitialEpoch + S{i}(2)/86400);
    AuxOrbit = AuxOrbit.Normalize(true, r0); 
    AuxOrbit = AuxOrbit.ChangeStateFormat('COE');
    AuxOrbit = AuxOrbit.DefineJ2Problem(J2, Re);
    AuxOrbit = AuxOrbit.AddPropagator('Osculating J2', Step);

    Constellation_1 = Constellation_1.AddOrbit(AuxOrbit);
end

%% Sensor network 
% Define an inertial observer
InitialState = [0 1 0];
InitialEpoch = juliandate(datetime('now'));
Sigma = diag([1e3 1e3 1e3]);
PD = 0.98;

Pc = 0.0;                   % Probability of false measurements
Vc = 10;                    % Number of false measurements per sensor (surveillance region)

InObs = Sensors.GibbsSensor(InitialEpoch, InitialState, Sigma, PD);

% Define a radar topocentric observer located at Madrid
InitialState = [deg2rad(0) deg2rad(-30)];
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
TelescopeObs = Sensors.TopocentricSensor(InitialEpoch, InitialState, Sigma, PD);

%% Observation process 
% Prepare the measurements
FinalObserveEpoch = EndEpoch;

InTime = [];
meas = [];
InState = [];

RadarTime = [];
meas_radar = [];
RadarState = [];

TelescopeTime = [];
meas_radec = [];
TelescopeState = [];

for i = 1:1
    [InTime_aux, meas_aux, InState_aux] = InObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);   

    InTime = [InTime; InTime_aux];
    meas = [meas; [i*ones(size(meas_aux,1),1) meas_aux]];
    InState = [InState; InState_aux];

    [RadarTime_aux, meas_radar_aux, RadarState_aux] = RadarObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

    RadarTime = [RadarTime; RadarTime_aux];
    meas_radar = [meas_radar; [i*ones(size(meas_radar_aux,1),1) meas_radar_aux]];
    RadarState = [RadarState; RadarState_aux];

    [TelescopeTime_aux, meas_radec_aux, TelescopeState_aux] = TelescopeObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

    TelescopeTime = [TelescopeTime; TelescopeTime_aux];
    meas_radec = [meas_radec; [i*ones(size(meas_radec_aux,1),1) meas_radec_aux]];
    TelescopeState = [TelescopeState; TelescopeState_aux];
end

ObservationSpan = [InTime; RadarTime; TelescopeTime];
[ObservationSpan, index] = sort(ObservationSpan);

Measurements = cell(length(ObservationSpan), 5);
Measurements(:,1) = num2cell(ObservationSpan);

for i = 1:length(index)
    if (index(i) <= size(InTime,1))
        Sigma = diag([1e4 1e4 1e4].^2/Re^2);
        Measurements(i,2) = { meas(index(i),:) };
        Measurements(i,3) = { InState(index(i),:) };
        Measurements(i,4) = { @(y)InObs.LikelihoodFunction(Sigma, meas(index(i),2:end).'/Re, y) };
        Measurements(i,5) = { @(y)InObs.ObservationProcess(InTime(index(i)), y, InState(index(i),:)) };

    elseif (index(i) <= size(InTime,1) + size(RadarTime,1))
        L = size(InTime,1);
        Sigma = diag([1e4 1e2].^2./[Re Re/Tc].^2);
        Measurements(i,2) = { meas_radar(index(i)-L,:) };
        Measurements(i,3) = { RadarState(index(i)-L,:) };
        Measurements(i,4) = { @(y)RadarObs.LikelihoodFunction(Sigma, meas_radar(index(i)-L,2:end).'./[Re; Re/Tc], y) };
        Measurements(i,5) = { @(y)RadarObs.ObservationProcess(RadarTime(index(i)-L), y, RadarState(index(i)-L,:)) };

    else
        Sigma = diag([deg2rad(1) deg2rad(1) deg2rad(0.1) deg2rad(0.1)].^2);
        L = size(InTime,1) + size(RadarTime,1);
        Measurements(i,2) = { meas_radec(index(i)-L,:) };
        Measurements(i,3) = { TelescopeState(index(i)-L,:) };
        Measurements(i,4) = { @(y)TelescopeObs.LikelihoodFunction(Sigma, meas_radec(index(i)-L,2:end).', y) };
        Measurements(i,5) = { @(y)TelescopeObs.ObservationProcess(TelescopeTime(index(i)-L), y, TelescopeState(index(i)-L,:)) };
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

%% True state
Re = 6378e3;
D(5) = sqrt(ElementSet(1)/Re); 
D(6) = sqrt((1-ElementSet(2)^2))* D(5);
D(7) = D(6) * cos(ElementSet(4));
D(1) = sin(ElementSet(4)/2) * cos((ElementSet(3)-ElementSet(5))/2);
D(2) = sin(ElementSet(4)/2) * sin((ElementSet(3)-ElementSet(5))/2);
D(3) = cos(ElementSet(4)/2) * sin((ElementSet(3)+ElementSet(5))/2);
D(4) = cos(ElementSet(4)/2) * cos((ElementSet(3)+ElementSet(5))/2);

%% Estimation: IOD
% Estimator configuration
IOD_filter = Filters.IOD_filter(30, 30, 5, PD, 1);

% Estimation
tic
[f, x, N_hat] = IOD_filter.BayesRecursion(ObservationSpan, Measurements);
running_time = toc;

%% Analysis
% Expectance of the number of targets in time 
N_error = N_hat - N;
Analysis.ExpectanceTargets = [mean(N_error) std(N_error)];

for i = 1:length(tspan)
    % CPEP
    Analysis.CPEP(i) = Filters.CPEP(N, N_hat, X, X_hat);

    % Hausdorff
    Analysis.HaussdorfDistance(i) = Filters.Hausdorff(X(:,i), X_hat(:,i));
end

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

% Constellation plot 
figure(10)
hold on
Constellation_1.OrbitSet{1,2}.PlotTrajectory(figure(10), InitialEpoch, EndEpoch);
for i = 1:size(Constellation_1.OrbitSet,1)
    Raux = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('ECI').Normalize(false,r0).StateEvolution(end,1:4);
    if (Raux(1) >= tspan(end))
        Raux = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('ECI').StateEvolution(end,2:4);
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