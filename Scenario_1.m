%% Combinatorial constellation macro-determination 
% Date: 02/02/2023
% Author: Sergio Cuevas del Valle

%% Constellation orbit determination. Scenario I

close all 
clear 

%% General user defined input
% Constants 
r0 = 6900e3;                % Characteristic distance of the Earth orbit    
mu = 3.86e14;               % Gravitional parameter of the Earth

Nmax = 4;                   % Number of targets
T = 800;                   % Number of orbital periods
tspan = 0:5e-1:T;             % Observation time span

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
Constellation_1 = Constellation('User defined', 1, 1, 1);
Constellation_1 = Constellation_1.ChangeTimeStep(1);

% Orbit definition
ElementType = 'COE'; 
ElementSet = [r0 1e-3 0 deg2rad(90) deg2rad(0)]; 

Sigma = 1e-5*eye(length(ElementSet));

for i = 1:size(S,1)
    % Generate a random anomaly 
    ElementSet = [ElementSet(1:5) 2*pi*rand()];

    % Add the orbit to the constellation
    AuxOrbit = Orbit(mu, ElementType, ElementSet, S{i}(1));
    AuxOrbit = AuxOrbit.SetFinalEpoch(S{i}(2));
    AuxOrbit = AuxOrbit.Normalize(true, r0); 
    AuxOrbit = AuxOrbit.ChangeStateFormat('COE');
    AuxOrbit = AuxOrbit.AddPropagator('Keplerian', 0.5);

    Constellation_1 = Constellation_1.AddOrbit(AuxOrbit);
end

% Set graphics
AuxOrbit.set_graphics();

% Compute the constellation parameters: number of planes, spacecraft and spacecraft per plane
Constellation_1.N = Constellation_1.NumberOfSpacecraft();
[Constellation_1.Np, Constellation_1.n] = Constellation_1.NumberOfPlanes();

% Propagation 
Constellation_1 = Constellation_1.Propagate(tspan(end));

%% Observation process 
% Measurements noise 
sigma_m = 1e-5;

% Random sampling of the state vector
s_meas = [];
time_meas = [];

for i = 1:size(Constellation_1.OrbitSet,1)
    AuxOrbit = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('COE').Normalize(false, r0);
    s_meas = [s_meas; AuxOrbit.StateEvolution(:,7)];
    time_meas = [time_meas; AuxOrbit.StateEvolution(:,1)];
end

s_meas = reshape(s_meas.',[size(s_meas,2)*size(s_meas,1) 1]); 
[time_meas, index] = sort(time_meas); 
s_meas = s_meas(index,:);

index = logical(randsrc(length(time_meas),1,[0, 1; 1-PD, PD]));
time_meas = time_meas(index);
s_meas = s_meas(index,:);

% Measurements
x = normrnd(cos(s_meas),sigma_m,size(s_meas,1),1);
y = normrnd(sin(s_meas),sigma_m,size(s_meas,1),1);

R = [x y]; 
meas = [time_meas R];

% Clutter generation 
index = logical(randsrc(length(tspan), 1, [0, 1; 1-Pc, Pc]));
theta_c = -2*pi+4*pi*rand(length(tspan),1);
x = normrnd(cos(theta_c),sigma_m,size(theta_c,1),1);
y = normrnd(sin(theta_c),sigma_m,size(theta_c,1),1);
clutter = [x y];
clutter = [tspan(index) clutter(index,:)];

if (size(clutter,1) > Vc)
    index = randi([0 size(clutter,1)], Vc,1); 
    clutter = clutter(index,:);
end

% Complete measurement set 
meas = [meas; clutter];
[~, index] = sort(meas(:,1)); 
meas = meas(index,:);

%% Estimator configuration
% UKF
UKF_estimator = EstimatorUKF('UKF-A', 2, 1E-2, 0); 
UKF_estimator = UKF_estimator.AssignStateProcess(1, @(time_step, theta)kinematic_proposal(mu, ElementSet(1), 0, time_step, theta, 0));
UKF_estimator = UKF_estimator.AssignObservationProcess(2, @(theta)radar(theta));
UKF_estimator = UKF_estimator.AdditiveCovariances(Sigma(1,1), sigma_m*eye(size(R,2)));
UKF_estimator = UKF_estimator.Init();

% EKF
EKF_estimator = EstimatorEKF;
EKF_estimator = EKF_estimator.AssignStateProcess(1, @(Q, time_step, theta, sigma)kinematic_proposal(mu, ElementSet(1), Q, time_step, theta, sigma));
EKF_estimator = EKF_estimator.AssignObservationProcess(2, @(theta)radar(theta));
EKF_estimator = EKF_estimator.AdditiveCovariances(Sigma(1,1), sigma_m*eye(size(R,2)));

% Mixture definition
J = 100;                            % Initial number of Gaussian components                                   
Jmax = 100;                         % Maximum number of Gaussian components
Mean = linspace(0,2*pi,J).';        % Random initialization of the mean of the components
Sigma = 1e-4*ones(J,1);                % Random initialization of the SD of the components

PHD = PHDFilter(J, Jmax, Mean, Sigma, PS, PD);

PHD = PHD.AddGaussian('Wrapped');           % Gaussian model to be used (in S^1)

PHD = PHD.DefineDomain(0,2*pi,1000);        % Domain definition and discretization
PHD = PHD.birth(1, pi, (0.5*pi)^2);         % Birth process definition
PHD = PHD.DefinePruning(1e-5, 4);           % Pruning process definition

% Likelihood process
PHD = PHD.AssignLikelihood(@(y,z,P)likelihood_function(y,z,P));

%% Estimation 
% Estimation
[f, X, N(2,:)] = PHD.kinematic_estimator(tspan, meas, EKF_estimator);
theta = PHD.Domain;

% State estimation 
if (~isempty(X{end}))
    Se = sort(X{end}(1,:));
else
    Se = [];
end

St = [];
for i = 1:size(Constellation_1.OrbitSet,1)
    AuxOrbit = Constellation_1.OrbitSet{i,2}.ChangeStateFormat('COE').Normalize(false, r0);
    St = [St; AuxOrbit.StateEvolution(end,end-1)];
end
St = mod(St,2*pi);
St = sort(St.');

%% Results
r = 1; 

figure
hold on
plot(theta, f(:,end), 'r')
for i = 1:size(St,2)
    xline(mod(St(i), 2*pi),'k'); 
end
xlabel('$\theta$')
ylabel('$p(\theta)$')
grid on;

x = r.*cos(theta);
y = r.*sin(theta);

figure
view(3)
hold on
plot3(r*cos(theta), r.*sin(theta), zeros(1,length(theta)))
plot3(x,y,f(:,end))
stem3(r*cos(St.'), r*sin(St.'), 10*ones(size(St)), 'k');
stem3(r.*cos(Se), r.*sin(Se), 10*ones(size(Se)), '*')
zlabel('$p(\theta)$')
grid on;

figure 
hold on
plot(tspan, N(1,:))
scatter(tspan, N(2,:))
xlabel('$t$')
ylabel('$N$')
legend('$\hat{N}$', '$N$')
grid on; 

if (0)
    figure
    view(3)
    [L, T] = meshgrid(t,theta);
    surface(L,T,f)
    xlabel('$t$')
    ylabel('$\theta$')
    zlabel('$f$')
    grid on;
end

if (0)
    dh = 50;
    W = figure;
    set(W, 'color', 'white');
    filename = 'Estimation.gif';
    view([-50 30]) 
    hold on

    xlabel('$x$');
    ylabel('$y$');
    zlabel('$p(\theta)$')
    grid on;

    % Mean anomaly manifold
    plot3(ElementSet(1).*cos(theta), ElementSet(1).*sin(theta), zeros(1,length(theta)));

    for i = 1:floor(length(tspan)/90):length(tspan)
        % Intensity function
        H = plot3(x, y, f(:,i), 'r');

        % States
        se = sort(X{i}(1,:));
        st = [];
        for j = 1:size(S,1)
            if (any(S{j}(:,1) == tspan(i)))
                index = S{j}(:,1) == tspan(i);
                pos = find(index, 1, 'first');
                st = [st; S{j}(pos,2)];
            end
        end

        J = scatter3(r*cos(st.'), r*sin(st.'), 10*ones(size(st)), 'k');
        V = scatter3(r.*cos(se), r.*sin(se), 10*ones(size(se)), '*');

        legend('$f(\theta)$', '$\mathbf{S}$', '$\hat{\mathbf{S}}$', 'AutoUpdate', 'off');

        drawnow;
        frame = getframe(W);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        if (i == 1) 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1e-2); 
        else 
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1e-2); 
        end 

        delete(J)
        delete(V)
        delete(H)
    end
    hold off
end

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