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
if (1)
    % Constants 
    r0 = 6928e3;                % Characteristic distance of the Earth orbit    
    mu = 3.986e14;              % Gravitional parameter of the Earth
    Re = 6378e3;                % Reference Earth radius
    Tc = sqrt(Re^2/mu);         % Characteristic time
    J2 = 1.08263e-3;            % Earth's J2 parameter
    
    % Constellation lifetime
    InitialEpoch = juliandate(datetime('now'));         % Initial epoch in JD
    T = 5;                                              % Number of days 
    EndEpoch = juliandate(datetime('now')+days(T));     % End epoch
    Step = 12 * 3600;                                    % Integration step in seconds
    tspan = 0:Step:T * 86400;                           % Relative lifetime in seconds
    
    % Target birth 
    PS = 1;                  % Probability of surviving
    PB = 0.000;              % Birth rate 
    
    % Constellation definition
    ElementType = 'COE'; 
    Plane(1,:) = Astrodynamics.sso_elements(10.5, [Re+600e3 1e-3 deg2rad(0) deg2rad(97) deg2rad(0)]);
    Plane(2,:) = Astrodynamics.sso_elements(13, [Re+400e3 1e-3 deg2rad(0) deg2rad(97) deg2rad(0)]);
    Plane(3,:) = [Re+500e3 1e-3 deg2rad(0) deg2rad(56) deg2rad(0) 0];
    Plane = [Plane; repmat(Plane(1,:), 7, 1); repmat(Plane(2,:), 6, 1); repmat(Plane(3,:), 5, 1); ];
    Nmax = 21;
    
    %% Target births and deaths 
    % Preallocation 
    N = size(Plane,1).*ones(1,length(tspan));    % Number of total spacecraft in time
    S = cell(size(Plane,1),1);                   % Time span of each target
    
    % Original births
    for i = 1:size(Plane,1)
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
        if (all(N >= size(Plane,1)))
            GoOn = false;
        end
    end
    
    %% Constellation definition (only one plane) 
    % Constellation definition 
    Constellation_1 = Constellation('User defined', size(Plane,1), 1, size(Plane,1));
    Constellation_1 = Constellation_1.ChangeTimeStep(Step);
    
    for i = 1:size(Plane,1)
        % Generate a random anomaly 
        if (i <= 8)
            ElementSet = [Plane(i,1:5) deg2rad(45) * (i-1)];
        elseif (i <= 15)
            ElementSet = [Plane(i,1:5) deg2rad(51.43) * (i-9)];
        else
            ElementSet = [Plane(i,1:5) deg2rad(60) * (i-16)];
        end
    
        % Add the orbit to the constellation
        AuxOrbit = Orbit(mu, ElementType, ElementSet, InitialEpoch + S{i}(1)/86400);
        AuxOrbit = AuxOrbit.SetFinalEpoch(InitialEpoch + S{i}(2)/86400);
        AuxOrbit = AuxOrbit.Normalize(true, Re); 
        AuxOrbit = AuxOrbit.ChangeStateFormat('COE');
        AuxOrbit = AuxOrbit.DefineJ2Problem(J2, Re);
        AuxOrbit = AuxOrbit.AddPropagator('Osculating J2', Step);
    
        Constellation_1 = Constellation_1.AddOrbit(AuxOrbit);
    end
    
    %% Sensor network 
    Pc = 0.0;                  % Probability of false measurements
    Vc = 0;                    % Number of false measurements per sensor (surveillance region)
    
    PD = 0.98;                 % Detection probability
    
    % Define an inertial observer
    InitialState = [0 1 0];
    InitialEpoch = juliandate(datetime('now'));
    Sigma = diag([5e1 5e1 5e1]);
    InObs = Sensors.GibbsSensor(InitialEpoch, InitialState, Sigma, PD);
    
    % Define a radar topocentric observer located at the Equator
    InitialState = [36.533 -6.283];
    InitialEpoch = juliandate(datetime('now'));
    Sigma = diag([70 10]);
    min_el = deg2rad(7);
    RadarObs = Sensors.RadarSensor(InitialEpoch, InitialState, Sigma, PD, min_el);
        
    % Define a telescope topocentric observer located at Canarias
    InitialState = [28.3 -15.59];
    InitialEpoch = juliandate(datetime('now'));
    Sigma = diag([deg2rad(0.001) deg2rad(0.001)]);
    min_el = deg2rad(7);
    TelescopeObs = Sensors.TopocentricSensor(InitialEpoch, InitialState, Sigma, PD, min_el, 'RADEC');
    
    %% Observation process 
    % Prepare the measurements
    FinalObserveEpoch = EndEpoch;
    
    InTime = [];
    InMeas = [];
    InState = [];
    
    RadarTime = [];
    RdMeas = [];
    RadarState = [];
    
    TelescopeTime = [];
    RaMeas = [];
    TelescopeState = [];
    
    for i = 1:size(Constellation_1.OrbitSet,1)
        % Inertial observer
        [time_aux, meas_aux, state_aux] = InObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);   
    
        InTime = [InTime; time_aux];
        InMeas = [InMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
        InState = [InState; state_aux];
        
        % Radar observer
        [time_aux, meas_aux, state_aux] = RadarObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);
    
        RadarTime = [RadarTime; time_aux];
        RdMeas = [RdMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
        RadarState = [RadarState; state_aux];
    
        % Telescope observer
%         [time_aux, meas_aux, state_aux] = TelescopeObs.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);
%     
%         TelescopeTime = [TelescopeTime; time_aux];
%         RaMeas = [RaMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
%         TelescopeState = [TelescopeState; state_aux];    
    end

    %% Pre-processing of the measurements
    ObservationSpan = [InTime; RadarTime; TelescopeTime];
    [ObservationSpan, index] = sort(ObservationSpan);
    
    Measurements = cell(length(ObservationSpan), 6);
    Measurements(:,1) = num2cell(ObservationSpan);
    
    for i = 1:length(index)
        if (index(i) <= size(InTime,1))
            Sigma = diag([1e3 1e3 1e3]./[Re Re Re]);

            L = 0;
            Measurements(i,2) = { InMeas(index(i)-L,:)./[1 Re Re Re] };
            Measurements(i,3) = { InState(index(i)-L,:) };
            Measurements(i,4) = { @(y)InObs.LikelihoodFunction(Sigma, InMeas(index(i)-L,2:end).'./[Re; Re; Re], y) };
            Measurements(i,5) = { @(y)InObs.ObservationProcess(InTime(index(i)-L), y, InState(index(i)-L,:)) };
            Measurements(i,6) = {reshape(Sigma, [], 1)};
            Measurements(i,7) = {'INERTIAL'};

        elseif (index(i) <= size(RadarTime,1) + size(InTime,1))
            Sigma = diag([5e2 5e1].^2 ./ [Re Re/Tc].^2);

            L = size(InTime,1);
            Measurements(i,2) = { RdMeas(index(i)-L,:) };
            Measurements(i,3) = { RadarState(index(i)-L,:) };
            Measurements(i,4) = { @(y)RadarObs.LikelihoodFunction(Sigma, RdMeas(index(i)-L,2:end).', y) };
            Measurements(i,5) = { @(y)RadarObs.ObservationProcess(RadarTime(index(i)-L), y, RadarState(index(i)-L,:)) };
            Measurements(i,6) = { reshape(Sigma, [], 1) };

            Measurements(i,7) = { 'RADAR' };

        else 
            Sigma = diag([deg2rad(0.001) deg2rad(0.001)]);

            L = size(RadarTime,1) + size(InTime,1);
            Measurements(i,2) = { RaMeas(index(i)-L,:) };
            Measurements(i,3) = { TelescopeState(index(i)-L,:) };
            Measurements(i,4) = { @(y)TelescopeObs.LikelihoodFunction(Sigma, RaMeas(index(i)-L,2:end).', y) };
            Measurements(i,5) = { @(y)TelescopeObs.ObservationProcess(TelescopeTime(index(i)-L), y, TelescopeState(index(i)-L,:)) };
            Measurements(i,6) = {reshape(Sigma, [], 1)};
            
            Measurements(i,7) = {'Telescope'};
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
        aux(k:length(ObservationSpan)) = repmat(aux(k-1), 1, length(k:length(ObservationSpan)) );
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

    save POC_IOD.mat;
else
    load POC_IOD.mat;

    %% Estimation: IOD
    % Estimator configuration
    PD = 0.98;
    UIOD_filter = Filters.UIOD_filter(1, 5e2, PS, PD);
    
    UIOD_filter.Lmin = 1.03; 
    UIOD_filter.Lmax = 2;
    UIOD_filter.emax = 0.2;

    % Estimation
    [X, N_hat, Prior, E_IOD] = UIOD_filter.BayesRecursion(ObservationSpan, Measurements);
    
    % IOD plane state
    D_hat = zeros(7,size(X{end}(1:7,:),2)); 
    coe_hat = zeros(size(X{end}(1:7,:),2), 6);
    for i = 1:size(X{end}(1:7,:),2)
        D_hat(:,i) = [X{end}(1:4,i); X{end}(5:7,i)];
        dumb = Astrodynamics.Delaunay2MyElements(D_hat(:,i), false);
        coe_hat(i,:) = Astrodynamics.Delaunay2COE(1, dumb, true).';
    end
    
    %% True plane state
    % Preallocation 
    D = cell(1,1);
    Dq = zeros(8, size(Constellation_1.OrbitSet,1));
    
    % Computation of the final 
    for i = 1:3
        % Delaunay sets
        coe = Constellation_1.OrbitSet{i,2}.StateEvolution(end,2:7);
        D{1}(:,i) = Astrodynamics.Delaunay2COE(1, coe, false).'; 
        
        % Quaternionic sets 
        Dq(:,i) = Astrodynamics.Delaunay2MyElements(D{1}(:,i), true);
    end
    
    MeanOrbit = Orbit(mu, ElementType, ElementSet, InitialEpoch + S{1}(1)/86400);
    MeanOrbit = MeanOrbit.SetFinalEpoch(InitialEpoch + S{1}(2)/86400);
    MeanOrbit = MeanOrbit.Normalize(true, Re); 
    MeanOrbit = MeanOrbit.ChangeStateFormat('COE');
    MeanOrbit = MeanOrbit.DefineJ2Problem(J2, Re);
    MeanOrbit = MeanOrbit.AddPropagator('Mean J2', Step);
    MeanOrbit = MeanOrbit.SetCurrentEpoch(InitialEpoch + S{1}(2)/86400);
    MeanOrbit = MeanOrbit.Propagate();
       
    %% Estimation: MTT-A
    % Extended tracking configuration 
    Gamma = 1;                  % Measurements per scan 
    
    % Estimator configuration
    iod_plane = [D_hat(:,1); 1e-5 * reshape(eye(6), [], 1)];
    MTT = Filters.MTT_filter(iod_plane, 1, 5e2, PS, PD, Gamma);
    MTT.ExtendedTarget = false;
    
    % Physical parameters
    MTT.mu = mu;
    MTT.epsilon = -J2; 
    MTT.Re = Re;
    
    % Anomaly partition
    MTT.nu = linspace(0,2*pi,1e3);
    
    MTT.Lmin = 1.03; 
    MTT.Lmax = 6;
    MTT.emax = 0.2;

    % Estimation
    [f, x, M_hat, ~, E] = MTT.BayesRecursion(ObservationSpan, Measurements);
    
end

%% Analysis
% Expectance of the number of targets in time 
N_error = N_hat-N;
Analysis.ExpectancePlanes = [mean(N_error) std(N_error)];

for i = 1:length(N_hat)
    % CPEP
    Analysis.CPEP(1,i) = Filters.Valuation.CPEP(Constellation_1.N, N_hat(i), X{i}(1:7,:), Dq(1:end-1,:), 0.1, 1);
    Analysis.CPEP(2,i) = Filters.Valuation.CPEP(Constellation_1.N, N_hat(i), X{i}(1:7,:), Dq(1:end-1,:), 0.1, 0);

    % Hausdorff
    Analysis.HaussdorfDistance(1,i) = Filters.Valuation.Hausdorff(X{i}(1:7,:), Dq(1:end-1,:), 1);
    Analysis.HaussdorfDistance(2,i) = Filters.Valuation.Hausdorff(X{i}(1:7,:), Dq(1:end-1,:), 0);
end

% Expectance of the number of targets in time 
M_error = M_hat-N;
Analysis.ExpectanceTargets = [mean(M_error) std(M_error)];

for i = 1:length(M_hat)
    % CPEP
    if (M_hat(i))
        Analysis.CPEP(3,i) = Filters.Valuation.CPEP(Constellation_1.N, M_hat(i), x{i}(1:7,:), Dq(1:end-1,:), 0.1, 1);
        Analysis.CPEP(4,i) = Filters.Valuation.CPEP(Constellation_1.N, M_hat(i), x{i}(1:7,:), Dq(1:end-1,:), 0.1, 0);
        Analysis.CPEP(5,i) = Filters.Valuation.CPEP(Constellation_1.n, M_hat(i), x{i}(8,:), Dq(end,:), 0.1, 2);
    
        % Hausdorff
        Analysis.HaussdorfDistance(3,i) = Filters.Valuation.Hausdorff(x{i}(1:7,:), Dq(1:end-1,:), 1);
        Analysis.HaussdorfDistance(4,i) = Filters.Valuation.Hausdorff(x{i}(1:7,:), Dq(1:end-1,:), 0);
        Analysis.HaussdorfDistance(5,i) = Filters.Valuation.Hausdorff(x{i}(8,:), Dq(end,:), 2);
    else
    end
end

%% Results
tspan = 86400 * ( ObservationSpan-ObservationSpan(1) );

figure 
scatter3(InMeas(:,2)/Re, InMeas(:,3)/Re, InMeas(:,4)/Re, 'filled');
grid on; 
xlabel('$X$')
ylabel('$Y$')
zlabel('$Z$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

figure 
scatter(RdMeas(:,2)/Re, RdMeas(:,3)*Tc/Re, 'filled');
grid on; 
xlabel('$\rho$')
ylabel('$\dot{\rho}$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure 
scatter(RaMeas(:,2), RaMeas(:,3), 'filled');
grid on; 
xlabel('$\alpha$')
ylabel('$\delta$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

%%
figure 
view(3)
hold on

no = [eye(3); zeros(1,3)];
nref = [eye(3); zeros(1,3)];
for j = 1:size(no,2)
    aux = QuaternionAlgebra.right_isoclinic(no(:,j)) * QuaternionAlgebra.quaternion_inverse(Dq(1:4,1));
    nref(:,j) = QuaternionAlgebra.right_isoclinic(Dq(1:4,1)) * aux;
end

quiver3(zeros(1,size(nref,2)), zeros(1,size(nref,2)), zeros(1,size(nref,2)), nref(1,:), nref(2,:), nref(3,:), 'b', 'LineWidth', 1.0); 


selected = [4];%[4 5 6 1 3];
n = no;
for i = 1:length(selected)
    for j = 1:size(no,2)
        aux = QuaternionAlgebra.right_isoclinic(no(:,j)) * QuaternionAlgebra.quaternion_inverse(D_hat(1:4,selected(i)));
        n(:,size(no,2)*(i-1)+j) = QuaternionAlgebra.right_isoclinic(D_hat(1:4,selected(i))) * aux;
    end
end
quiver3(zeros(1,size(n,2)), zeros(1,size(n,2)), zeros(1,size(n,2)), n(1,:), n(2,:), n(3,:), 'r', 'LineWidth', 1.0); 

legend('$\mathcal{P}$', '$\hat{\mathcal{P}}$', 'AutoUpdate', 'off', 'Location', 'best')

m = 100;
[aa, bb, cc] = sphere(m);
h = surf(aa, bb, cc);
set(h, 'FaceColor', [0 0 1], 'FaceAlpha', 0.00001, 'FaceLighting', 'gouraud', 'EdgeColor', 'none')

hold off
xlabel('$X$')
ylabel('$Y$')
zlabel('$Z$')
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

%%
figure
hold on
yline(3, 'b')
scatter(tspan/3600, N_hat, 'r', 'filled')
hold off
legend('$N_p$',' $\hat{N}_p$');
xlabel('Epoch $t$ [h]')
ylabel('N. planes')
grid on;

%%
figure
hold on
plot(tspan, E_IOD, 'o-')
hold off
xlabel('Epoch $t$')
ylabel('Diff. entropy $E_{max}$')
grid on;

%%
figure
hold on
plot(tspan/3600, Analysis.CPEP(1,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{CPEP}_{\mathbf{q}}$')
grid on;

figure
hold on
plot(tspan/3600, Analysis.CPEP(2,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{CPEP}_{\mathbf{\Gamma}}$')
grid on;

figure
hold on
plot(tspan/3600, Analysis.HaussdorfDistance(1,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{H}_{\mathbf{q}}$')
grid on;

figure
hold on
plot(tspan/3600, Analysis.HaussdorfDistance(2,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{H}_{\mathbf{\Gamma}}$')
grid on;

%%
figure
hold on
plot(tspan/3600, Analysis.CPEP(3,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{CPEP}_{\mathbf{q}}$')
grid on;

figure
hold on
plot(tspan/3600, Analysis.CPEP(4,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{CPEP}_{\mathbf{\Gamma}}$')
grid on;

figure
hold on
plot(tspan/3600, Analysis.CPEP(5,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{CPEP}_{M}$')
grid on;

figure
hold on
plot(tspan/3600, Analysis.HaussdorfDistance(3,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{H}_{\mathbf{q}}$')
grid on;

figure
hold on
plot(tspan/3600, Analysis.HaussdorfDistance(4,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{H}_{\mathbf{\Gamma}}$')
grid on;

figure
hold on
plot(tspan/3600, Analysis.HaussdorfDistance(5,:), 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('$\textrm{H}_{M}$')
grid on;

%% 
% Anomaly plot 
pos = length(M_hat)-1;
figure 
view(3)
hold on
plot3(cos(MTT.nu), sin(MTT.nu), f{pos}); 
for i = 1:min(size(x{pos},2), 1)
    stem3(cos(x{pos}(8,i)), sin(x{pos}(8,i)), 1, 'r', 'filled');
end
xlabel('$X$')
ylabel('$Y$')
for i = 1:min(1,size(Constellation_1.OrbitSet,1))
    if (Constellation_1.OrbitSet{i,2}.InitialEpoch + S{i}(end)/86400 >= ObservationSpan(pos))
        diff_time = Constellation_1.OrbitSet{i,2}.InitialEpoch + Constellation_1.OrbitSet{i,2}.Normalize(false, r0).StateEvolution(:,1)/86400 >= ObservationSpan(pos);
        index = find(diff_time, 1, 'first');
        Set = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('COE').StateEvolution(index,7);
        scatter3(cos(Set), sin(Set), 1, 'filled', 'b');
    end
end
legend('$f(M)$', '$\hat{M}_i$', '$M$', 'AutoUpdate', 'off')
plot(cos(MTT.nu), sin(MTT.nu), 'k');
for i = 2:size(x{pos},2)
    stem3(cos(x{pos}(8,i)), sin(x{pos}(8,i)), 1, 'r', 'filled');
end
for i = 2:size(Constellation_1.OrbitSet,1)
    if (Constellation_1.OrbitSet{i,2}.InitialEpoch + S{i}(end)/86400 >= ObservationSpan(pos))
        diff_time = Constellation_1.OrbitSet{i,2}.InitialEpoch + Constellation_1.OrbitSet{i,2}.Normalize(false, r0).StateEvolution(:,1)/86400 >= ObservationSpan(pos);
        index = find(diff_time, 1, 'first');
        Set = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('COE').StateEvolution(index,7);
        scatter3(cos(Set), sin(Set), 1, 'filled', 'b');
    end
end
grid on; 
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

%%
figure
hold on
plot(tspan/3600, E, 'o-')
hold off
xlabel('Epoch $t$ [h]')
ylabel('Diff. entropy $E_{max}$')
grid on;

%%
figure
hold on
plot(tspan/3600, N)
scatter(tspan/3600, M_hat, 'filled')
hold off
legend('$n_p$',' $\hat{n}_p$');
xlabel('Epoch $t$ [h]')
ylabel('N. spacecraft')
grid on;
ylim([0 4])
