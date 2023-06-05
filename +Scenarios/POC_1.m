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
J2 = 1.08263e-3;            % Earth's J2 parameter
Nmax = 2;                   % Number of targets

% Constellation lifetime
InitialEpoch = juliandate(datetime('now'));         % Initial epoch in JD
T = 1;                                              % Number of days 
EndEpoch = juliandate(datetime('now')+days(T));     % End epoch
Step = 1200;                                         % Integration step in seconds
tspan = 0:Step:T * 86400;                           % Relative lifetime in seconds

% Target birth 
PS = 1;                  % Probability of surviving
PB = 0.000;                 % Birth rate 

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
Plane(1,:) = [r0 1e-2 deg2rad(20) deg2rad(15) deg2rad(10)];
Plane(2,:) = [r0 1e-3 deg2rad(20) deg2rad(15) deg2rad(10)]; 
Plane(3,:) = [r0 1e-3 deg2rad(20) deg2rad(15) deg2rad(10)];

for i = 1:size(S,1)
    % Generate a random anomaly 
    ElementSet = [Plane(i,1:5) 2*pi*rand()];

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
Pc = 0.0;                   % Probability of false measurements
Vc = 10;                    % Number of false measurements per sensor (surveillance region)

PD = 0.98;                  % Detection probability

% Define an inertial observer
InitialState = [0 1 0];
InitialEpoch = juliandate(datetime('now'));
Sigma = diag([1e1 1e1 1e1]);
InObs = Sensors.GibbsSensor(InitialEpoch, InitialState, Sigma, PD);

% Define an anomaly observer
InitialState = 0;
Sigma = deg2rad(0.001)^2;
AnObs = Sensors.AnomalySensor(InitialEpoch, InitialState, Sigma, PD);

% Define an Delaunay observer
InitialState = 0;
Sigma = deg2rad(0.001)^2 * eye(6);
DyObs = Sensors.DelaunaySensor(InitialEpoch, InitialState, Sigma, PD);

% Define a radar topocentric observer located at the Equator
InitialState = [0 -30];
InitialEpoch = juliandate(datetime('now'));
Sigma = eye(2,2);
FOV = deg2rad(40);
RadarObs = Sensors.RadarSensor(InitialEpoch, InitialState, Sigma, PD, FOV);

% Define a telescope topocentric observer located at Madrid
InitialState = [40 -3];
InitialEpoch = juliandate(datetime('now'));
Sigma = diag([deg2rad(0.01) deg2rad(0.01) deg2rad(0.001) deg2rad(0.001)]);
FOV = deg2rad(40);
TelescopeObs = Sensors.TopocentricSensor(InitialEpoch, InitialState, Sigma, PD, FOV);

% Define a telescope topocentric observer located at the north pole
InitialState = [0 90];
InitialEpoch = juliandate(datetime('now'));
FOV = deg2rad(40);
% Sigma = diag([deg2rad(0.01) deg2rad(0.01) deg2rad(0.001) deg2rad(0.001)]);
Sigma = diag([deg2rad(0.001) deg2rad(0.001)]);
TelescopeObs2 = Sensors.TopocentricSensor(InitialEpoch, InitialState, Sigma, PD, FOV, 'RADEC');

%% Observation process 
% Prepare the measurements
FinalObserveEpoch = EndEpoch;

InTime = [];
InMeas = [];
InState = [];

DyTime = [];
DyMeas = [];
DyState = [];

AnTime = [];
AnMeas = [];
AnState = [];

RadarTime = [];
RdMeas = [];
RadarState = [];

TelescopeTime = [];
RaMeas = [];
TelescopeState = [];

TelescopeTime2 = [];
RaMeas2 = [];
TelescopeState2 = [];

for i = 1:size(Constellation_1.OrbitSet,1)
    % Inertial observer
    [time_aux, meas_aux, state_aux] = InObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);   

    InTime = [InTime; time_aux];
    InMeas = [InMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
    InState = [InState; state_aux];

    % Delaunay observer
    [time_aux, meas_aux, state_aux] = DyObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);   

    DyTime = [DyTime; time_aux];
    DyMeas = [DyMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
    DyState = [DyState; state_aux];

    % Anomaly observer
    [time_aux, meas_aux, state_aux] = AnObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);   

    AnTime = [AnTime; time_aux];
    AnMeas = [AnMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
    AnState = [AnState; state_aux];

    % Radar observer
    [time_aux, meas_aux, state_aux] = RadarObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

    RadarTime = [RadarTime; time_aux];
    RdMeas = [RdMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
    RadarState = [RadarState; state_aux];

    % Telescope observer
    [time_aux, meas_aux, state_aux] = TelescopeObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

    TelescopeTime = [TelescopeTime; time_aux];
    RaMeas = [RaMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
    TelescopeState = [TelescopeState; state_aux];

    % Second telescope observer
    [time_aux, meas_aux, state_aux] = TelescopeObs2.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

    TelescopeTime2 = [TelescopeTime2; time_aux];
    RaMeas2 = [RaMeas2; [i*ones(size(meas_aux,1),1) meas_aux]];
    TelescopeState2 = [TelescopeState2; state_aux];
end

% ObservationSpan = [InTime; DyTime; AnTime; RadarTime; TelescopeTime; TelescopeTime2];
ObservationSpan = [DyTime];
[ObservationSpan, index] = sort(ObservationSpan);

Measurements = cell(length(ObservationSpan), 6);
Measurements(:,1) = num2cell(ObservationSpan);

for i = 1:length(index)
%     if (index(i) <= size(InTime,1))
%         Sigma = diag([1e5 1e5 1e5].^2/Re^2);
%         Measurements(i,2) = { meas(index(i),:)./[1 Re Re Re] };
%         Measurements(i,3) = { InState(index(i),:) };
%         Measurements(i,4) = { @(y)InObs.LikelihoodFunction(Sigma, meas(index(i),2:end).'/Re, y) };
%         Measurements(i,5) = { @(y)InObs.ObservationProcess(InTime(index(i)), y, InState(index(i),:)) };
%         Measurements(i,6) = { reshape(Sigma, [], 1) };
%         Measurements(i,7) = { 'INERTIAL' };
%     end
% 
%     elseif (index(i) <= size(InTime,1) + size(RadarTime,1))
%         L = size(InTime,1);
%         Sigma = diag([1e3 1e2].^2./[Re Re/Tc].^2);
%         Measurements(i,2) = { meas_radar(index(i)-L,:)./[1 Re Re/Tc] };
%         Measurements(i,3) = { RadarState(index(i)-L,:) };
%         Measurements(i,4) = { @(y)RadarObs.LikelihoodFunction(Sigma, meas_radar(index(i)-L,2:end).'./[Re; Re/Tc], y) };
%         Measurements(i,5) = { @(y)RadarObs.ObservationProcess(RadarTime(index(i)-L), y, RadarState(index(i)-L,:)) };
%         Measurements(i,6) = {reshape(Sigma, [], 1)};
%         Measurements(i,7) = {'RADAR'};
% 
%     elseif (index(i) <= size(InTime,1) + size(RadarTime,1) + size(TelescopeTime,1))
%         Sigma = diag([deg2rad(1) deg2rad(1) deg2rad(0.1) deg2rad(0.1)].^2);
%         L = size(InTime,1) + size(RadarTime,1);
%         Measurements(i,2) = { meas_radec(index(i)-L,:) };
%         Measurements(i,3) = { TelescopeState(index(i)-L,:) };
%         Measurements(i,4) = { @(y)TelescopeObs.LikelihoodFunction(Sigma, meas_radec(index(i)-L,2:end).', y) };
%         Measurements(i,5) = { @(y)TelescopeObs.ObservationProcess(TelescopeTime(index(i)-L), y, TelescopeState(index(i)-L,:)) };
%         Measurements(i,6) = {reshape(Sigma, [], 1)};
%         Measurements(i,7) = {'Telescope'};
% 
%     else
% 
%     end

    if (0)
        Sigma = diag([deg2rad(0.1) deg2rad(0.1) deg2rad(0.01) deg2rad(0.01)]);
        Sigma = diag([deg2rad(1e-1) deg2rad(1e-1)]);
        Measurements(i,2) = { RaMeas2(index(i),:) };
        Measurements(i,3) = { TelescopeState2(index(i),:) };
        Measurements(i,4) = { @(y)TelescopeObs2.LikelihoodFunction(Sigma, RaMeas2(index(i),2:end).', y) };
        Measurements(i,5) = { @(y)TelescopeObs2.ObservationProcess(TelescopeTime2(index(i)), y, TelescopeState2(index(i),:)) };
        Measurements(i,6) = {reshape(Sigma, [], 1)};
        Measurements(i,7) = {'Telescope'};
    end

    if (0) 
        Sigma = deg2rad(1e-1).^2;
        Measurements(i,2) = { AnMeas(index(i),:) };
        Measurements(i,3) = { AnState(index(i),:) };
        Measurements(i,4) = { @(y)AnObs.LikelihoodFunction(Sigma, AnMeas(index(i),2:end).', y) };
        Measurements(i,5) = { @(y)AnObs.ObservationProcess(AnTime(index(i)), y, AnState(index(i),:)) };
        Measurements(i,6) = { reshape(Sigma, [], 1) };
        Measurements(i,7) = { 'ANOMALY' };
    end

    if (1)
        Sigma = blkdiag(deg2rad(5)^2 * eye(3), 5e-2 * eye(3));
        Measurements(i,2) = { DyMeas(index(i),:) };
        Measurements(i,3) = { DyState(index(i),:) };
        Measurements(i,4) = { @(y)DyObs.LikelihoodFunction(Sigma, DyMeas(index(i),2:end).', y) };
        Measurements(i,5) = { @(y)DyObs.ObservationProcess(DyTime(index(i)), y, DyState(index(i),:)) };
        Measurements(i,6) = { reshape(Sigma, [], 1) };
        Measurements(i,7) = { 'DELAUNAY' };
    end

    if (0)
        Sigma = diag([1e2 1e2].^2 ./ [Re Re/Tc].^2);
        Measurements(i,2) = { RdMeas(index(i),:) ./ [1 Re Re/Tc] };
        Measurements(i,3) = { RadarState(index(i),:) };
        Measurements(i,4) = { @(y)RadarObs.LikelihoodFunction(Sigma, RdMeas(index(i),2:end).', y) };
        Measurements(i,5) = { @(y)RadarObs.ObservationProcess(RadarTime(index(i)), y, RadarState(index(i),:)) };
        Measurements(i,6) = { reshape(Sigma, [], 1) };
        Measurements(i,7) = { 'RADAR' };
    end
end

i = 1;
k = 1;
aux = zeros(1,length(ObservationSpan));
while (i < length(tspan) && k <= length(ObservationSpan))
    time_diff = tspan(i)/86400 + InitialEpoch <= ObservationSpan & tspan(i+1)/86400 + InitialEpoch > ObservationSpan;
    count = sum(time_diff);
    aux(k:k-1+count) = repmat(N(i), 1, count);
    i = i + 1;
    k = k + count;
end

if (k < length(ObservationSpan))
    aux(k:length(ObservationSpan)) = repmat(aux(k-1), 1, length(ObservationSpan)-k);
end

N = aux;

%% Propagation 
% Propagation 
Constellation_1 = Constellation_1.Propagate(EndEpoch);

% Set graphics
AuxOrbit.set_graphics();

% Compute the constellation parameters: number of planes, spacecraft and spacecraft per plane
Constellation_1.N = Constellation_1.NumberOfSpacecraft();
[Constellation_1.Np, Constellation_1.n] = Constellation_1.NumberOfPlanes();

%% True plane state
D = Astrodynamics.Delaunay2COE(1, [ElementSet(1:5) 0] ./ [Re 1 1 1 1 1], false).'; 

%% Estimation: IOD
if (0)
    % Estimator configuration
    UIOD_filter = Filters.UIOD_filter(1, 2e2, PS, PD);
    
    % Estimation
    [x_IOD, N_hat, Prior, E_IOD] = UIOD_filter.BayesRecursion(ObservationSpan, Measurements);
    
    % IOD plane state
    D_hat = Astrodynamics.Delaunay2MyElements(x_IOD{end}(1:7,1), false);
else
    D_hat = Astrodynamics.Delaunay2MyElements(D, true);
end

D_hat = [D_hat(1:end-1,1); 1e-7 * reshape(eye(6), [], 1)];

%% Estimation: tracking 
% Extended tracking configuration 
Gamma = 1;                  % Measurements per scan 

% Estimator configuration
if (1)
    MTT = Filters.MTT_filter(D, 1, 5e1, PS, PD, Gamma);
else
    MTT = Filters.CMTT_filter(D, 1, 1e2, PS, PD, Gamma);
end
MTT.ExtendedTarget = false;
MTT.planes = D_hat;

% Physical parameters
MTT.mu = mu;
MTT.epsilon = -J2; 
MTT.Re = Re;

% Anomaly partition
MTT.nu = linspace(0,2*pi,1e3);

% Estimation
tic
[f, x, N_hat, Prior, E] = MTT.BayesRecursion(ObservationSpan, Measurements);
running_time = toc;

%% Analysis
% Expectance of the number of targets in time 
N_error = N_hat-N;
Analysis.ExpectanceTargets = [mean(N_error) std(N_error)];

% for i = 1:length(tspan)
%     % CPEP
%     Analysis.CPEP(i) = Filters.CPEP(Constellation_1.n, N_hat, X, X_hat);
% 
%     % Hausdorff
%     Analysis.HaussdorfDistance(i) = Filters.Hausdorff(X(:,i), X_hat(:,i));
% end

%% Results
tspan = 86400 * ( ObservationSpan-ObservationSpan(1) );
pos = size(N_hat,2)-7;

% figure 
% view(3)
% no = [0;0;1];
% for i = 1:size(x_IOD{pos},2)
%     aux = QuaternionAlgebra.right_isoclinic([no; 0]) * QuaternionAlgebra.quaternion_inverse(x_IOD{pos}(1:4,i));
%     n(:,i) = QuaternionAlgebra.right_isoclinic(x_IOD{pos}(1:4,i)) * aux;
% end
% hold on
% m = 100;
% [aa, bb, cc] = sphere(m);
% h = surf(aa, bb, cc);
% set(h, 'FaceColor',[0 0 1], 'FaceAlpha',0.1,'FaceLighting','gouraud','EdgeColor','none')
% scatter3(n(1,:), n(2,:), n(3,:), 'filled'); 
% hold off
% xlabel('$X$')
% ylabel('$Y$')
% ylabel('$Z$')
% grid on;

figure
hold on
plot(tspan, N)
scatter(tspan, N_hat); 
hold off
legend('$N_p$',' $\hat{N}_p$');
xlabel('Epoch $t$')
ylabel('N. planes')
grid on;

figure
hold on
plot(tspan, N)
scatter(tspan, N_hat); 
hold off
legend('$N_p$',' $\hat{N}_p$');
xlabel('Epoch $t$')
ylabel('N. planes')
grid on;

figure
hold on
plot(tspan, E)
hold off
xlabel('Epoch $t$')
ylabel('Diff. entropy $E_{max}$')
grid on;

% Constellation plot 
figure(10)
hold on
Constellation_1.OrbitSet{1,2}.PlotTrajectory(figure(10), Constellation_1.OrbitSet{1,2}.InitialEpoch, Constellation_1.OrbitSet{1,2}.FinalEpoch);
for i = 1:size(Constellation_1.OrbitSet,1)
    Raux = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('ECI').Normalize(false,r0).StateEvolution(end,1);
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
figure 
view(3)
hold on
plot3(cos(MTT.nu), sin(MTT.nu), f{pos}); 
plot(cos(MTT.nu), sin(MTT.nu), 'k'); 
for i = 1:min(size(x{pos},2), 1)
    stem3(cos(x{pos}(8,i)), sin(x{pos}(8,i)), 1, 'k');
end
xlabel('$X$')
ylabel('$Y$')
for i = 1:min(1,size(Constellation_1.OrbitSet,1))
    if (Constellation_1.OrbitSet{i,2}.InitialEpoch + S{i}(end)/86400 >= ObservationSpan(pos))
        diff_time = Constellation_1.OrbitSet{i,2}.InitialEpoch + Constellation_1.OrbitSet{i,2}.Normalize(false, r0).StateEvolution(:,1)/86400 >= ObservationSpan(pos);
        index = find(diff_time, 1, 'first');
        Set = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('COE').StateEvolution(index,7);
        scatter3(cos(Set), sin(Set), 1, 'filled', 'r');
    end
end
legend('$f(M)$', '$\hat{M}_i$', '$M$', 'AutoUpdate', 'off')

for i = 2:size(x{pos},2)
    stem3(cos(x{pos}(8,i)), sin(x{pos}(8,i)), 1, 'k');
end
for i = 2:size(Constellation_1.OrbitSet,1)
    if (Constellation_1.OrbitSet{i,2}.InitialEpoch + S{i}(end)/86400 >= ObservationSpan(pos))
        diff_time = Constellation_1.OrbitSet{i,2}.InitialEpoch + Constellation_1.OrbitSet{i,2}.Normalize(false, r0).StateEvolution(:,1)/86400 >= ObservationSpan(pos);
        index = find(diff_time, 1, 'first');
        Set = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('COE').StateEvolution(index,7);
        scatter3(cos(Set), sin(Set), 1, 'filled', 'r');
    end
end
grid on; 
