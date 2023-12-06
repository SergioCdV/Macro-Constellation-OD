%% Combinatorial constellation macro-determination 
% Date: 0/02/2023
% Author: Sergio Cuevas del Valle

%% Constellation orbit determination. Scenario IOD %%
% This script provides the IOD scenario in which to test the on-manifold
% PHD to identify individual targets %

close all 
clear 
rng(1);

%% General user defined input
if (0)
    % Constants 
    r0 = 6928e3;                % Characteristic distance of the Earth orbit    
    mu = 3.986e14;              % Gravitional parameter of the Earth
    Re = 6378e3;                % Reference Earth radius
    Tc = sqrt(mu/Re^3);         % Characteristic time
    J2 = 1.08263e-3;            % Earth's J2 parameter
    
    % Constellation lifetime
    InitialEpoch = juliandate(datetime('now'));         % Initial epoch in JD
    T = 1;                                              % Number of days 
    EndEpoch = juliandate(datetime('now')+days(T));     % End epoch
    Step = 0.16 * 3600;                                 % Integration step in seconds
    tspan = 0:Step:T * 86400;                           % Relative lifetime in seconds
    
    % Target birth 
    PS = 1;                  % Probability of surviving
    PB = 0.000;              % Birth rate 
    
    % Constellation definition
    ElementType = 'COE'; 
    Plane(1,:) = Astrodynamics.sso_elements(10.5, [Re+600e3 1e-3 deg2rad(0) deg2rad(97) deg2rad(0)]);
    Plane(2,:) = Astrodynamics.sso_elements(13,   [Re+400e3 1e-3 deg2rad(0) deg2rad(97) deg2rad(0)]);
    Plane(3,:) = [Re+500e3 1e-3 deg2rad(0) deg2rad(56) deg2rad(0) 0];
    Plane = [Plane; repmat(Plane(1,:), 1, 1); repmat(Plane(2,:), 1, 1); repmat(Plane(3,:), 1, 1); ];
    Nmax = 3;
    
    %% Target births and deaths 
    % Preallocation 
    N = size(Plane,1) .* ones(1,length(tspan));    % Number of total spacecraft in time
    S = cell(size(Plane,1),1);                     % Time span of each target
    
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
    PD = 0.98;                 % Probability of detection of the spacecraft
    sampling_rate = 10000;     % Sampling rate [s]

    % Define Nr radar topocentric observer located at the Equator (evenly
    % distributed)
    nr = 12;                                                    % Number of radars
    lat = zeros(nr, 1);                                         % Latitude of the radars [deg] 
    long = linspace(-180, 180, nr+1).';                         % Longitude of the radars [deg]
    long = long(1:end-1);
    
    InitialState = [lat long];
    InitialEpoch = juliandate(datetime('now'));
    Sigma = diag([70 1]);
    min_el = deg2rad(3);

    for i = 1:nr
        RadarObs{i} = Sensors.RadarSensor(InitialEpoch, InitialState(i,:), Sigma, PD, min_el);
    end
        
    % Define a telescope topocentric observer located at Canarias
    InitialEpoch = juliandate(datetime('now'));
    Sigma = diag([deg2rad(0.001) deg2rad(0.001)]);
    min_el = deg2rad(7);
    
    % Define several telescope topocentric observers
    InitialState = [+28.3 -15.59; ...
                    +30.70 34.8; ...
                    +20.70 -156.3571; ...
                    +30.70 -104.02; ...
                    -32.37 20.81; ...
                    -2.89 53.189; ...
                    -30.167 -70.800; ...
                    -31.272 149.0616; ...
                    +32.2 80.0];

    for i = 1:size(InitialState,1)
        TelescopeObs{i} = Sensors.TopocentricSensor(InitialEpoch, InitialState, Sigma, PD, deg2rad(7), 'RADEC');
    end
    
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
        % Radar observer
        for j = 1:length(RadarObs)
            [time_aux, meas_aux, state_aux] = RadarObs{j}.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);

            % Sample the measurement set based on the sampling rate
%             diff_time = diff(time_aux);
%             time_aux = time_aux(diff_time > sampling_rate);
%             meas_aux = meas_aux(diff_time > sampling_rate,:);
%             state_aux = state_aux(diff_time > sampling_rate,:);

            RadarTime = [RadarTime; time_aux];
            RdMeas = [RdMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
            RadarState = [RadarState; state_aux];
        end
    
        % Telescope observer
        for j = 1:length(TelescopeObs)
            [time_aux, meas_aux, state_aux] = TelescopeObs{j}.Observe(Constellation_1.OrbitSet{i,2}, FinalObserveEpoch);
        
            % Sample the measurement set based on the sampling rate
%             diff_time = diff(time_aux);
%             time_aux = time_aux(diff_time > sampling_rate);
%             meas_aux = meas_aux(diff_time > sampling_rate,:);
%             state_aux = state_aux(diff_time > sampling_rate,:);

            TelescopeTime = [TelescopeTime; time_aux];
            RaMeas = [RaMeas; [i*ones(size(meas_aux,1),1) meas_aux]];
            TelescopeState = [TelescopeState; state_aux];    
        end
    end

    %% Pre-processing of the measurements
    ObservationSpan = [RadarTime; TelescopeTime];
    [ObservationSpan, index] = sort(ObservationSpan);
    
    Measurements = cell(length(ObservationSpan), 6);
    Measurements(:,1) = num2cell(ObservationSpan);
    
    for i = 1:length(index)
        if ( index(i) <= size(RadarTime,1) )
            Sigma = diag([5e2 5e1].^2 ./ [Re Re/Tc].^2);

            L = size(InTime,1);
            Measurements(i,2) = { RdMeas(index(i)-L,:) };
            Measurements(i,3) = { RadarState(index(i)-L,:) };
            Measurements(i,4) = { @(y)RadarObs{1}.LikelihoodFunction(Sigma, RdMeas(index(i)-L,2:end).', y) };
            Measurements(i,5) = { @(y)RadarObs{1}.ObservationProcess(RadarTime(index(i)-L), y, RadarState(index(i)-L,:)) };
            Measurements(i,6) = { reshape(Sigma, [], 1) };

            Measurements(i,7) = { 'RADAR' };

        else 
            Sigma = diag([deg2rad(0.001) deg2rad(0.001)]);

            L = size(RadarTime,1) + size(InTime,1);
            Measurements(i,2) = { RaMeas(index(i)-L,:) };
            Measurements(i,3) = { TelescopeState(index(i)-L,:) };
            Measurements(i,4) = { @(y)TelescopeObs{1}.LikelihoodFunction(Sigma, RaMeas(index(i)-L,2:end).', y) };
            Measurements(i,5) = { @(y)TelescopeObs{1}.ObservationProcess(TelescopeTime(index(i)-L), y, TelescopeState(index(i)-L,:)) };
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

    save SciTech2024_IOD.mat;
else
    load SciTech2024_IOD.mat;

    %% Estimation: IOD
    % Estimator configuration
    PD = 0.98;
    PHIOD_filter = Filters.PHIOD_filter(1, 1e3, PS, PD);
    
    PHIOD_filter.Lmin = 1.03; 
    PHIOD_filter.Lmax = 2;
    PHIOD_filter.emax = 0.2;

    % Estimation
    [X, N_hat, Prior, E_IOD] = PHIOD_filter.BayesRecursion(ObservationSpan, Measurements);
    
    % IOD spacecraft states
    D_hat = zeros(8,size(X{end}(1:8,:),2)); 
    coe_hat = zeros(size(X{end}(1:8,:),2), 6);
    for i = 1:size(X{end}(1:8,:),2)
        D_hat(:,i) = [X{end}(1:4,i); X{end}(5:8,i)];
        dumb = Astrodynamics.Delaunay2MyElements(D_hat(:,i), false);
        coe_hat(i,:) = Astrodynamics.Delaunay2COE(1, dumb, true).';
    end
    
    %% True plane state
    % Preallocation 
    D = cell(1,1);
    Dq = zeros(8, size(Constellation_1.OrbitSet,1));
    
    % Computation of the final element sets
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

%% Results
tspan = 86400 * ( ObservationSpan-ObservationSpan(1) );

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