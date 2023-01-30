%% Combinatorial constellation macro-determination 
% Date: 27/10/2022
% Author: Sergio Cuevas del Valle

%% Toy model for a multi-state Bayesian filter implementation
% The dynamics of a single spacecraft in a equatorial, planar circular
% orbit are analysed, based on inertial measurements of the radial vector

set_graphics();
close all; 
clear;

%% Input data 
% Constants
mu = 1;                 % Gravitational parameter

T = 10;                 % Stopping simulation time
tspan = 0:5e-1:2*pi*T; 

% Initial conditions
Nmax = 4;               % Number of constellation plane targets

% Distribution of the geometrical set (classical orbital elements)
sigma_oe = [0.5 1e-4 2*pi 2*pi 2*pi];
mu_oe = [2 1e-3 0 0 0];

%% Generation of targets
% Births 
PS = 0.999;                                                 % Probability of surviving
PB = 0.05;                                                  % Birth rate

% Target generation and death
N = Nmax.*ones(1,length(tspan));

S = cell(Nmax,1);
for i = 1:Nmax
    % Initial conditions
    theta = 2*pi*rand;
    set = normrnd(mu_oe, sigma_oe);

    % Deaths of the original plane 
    deaths = logical(randsrc(length(tspan),1,[0, 1; PS, 1-PS]));
    td = find(deaths, 1, 'first'); 

    % Final states
    if (isempty(td))
        td = length(tspan);
    else
        N(1,td:end) = N(1,td:end)-1;
    end
    S{i} = [0 tspan(td) set theta];           
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
        theta = 2*pi*rand;
        set = normrnd(mu_oe, sigma_oe);
        S{i+Nmax} = [tspan(tb) tspan(td) set theta];          % Targets 
        i = i+1;
        N(1,tb:td) = N(1,tb:td)+1; 
        if (max(N) > Nmax)
            a = 1;
        end
    end

    % Convergence 
    if (all(N >= Nmax))
        GoOn = false;
    end
end

%% True motion integration and measurements generation 
% True orbital dynamics
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);

for i = 1:size(S,1)
    [t, s] = ode45(@(t,s)dynamics(mu, t, s), tspan(tspan >= S{i}(1) & tspan <= S{i}(2)), S{i}(4:end), options);
    S{i} = [t s];
end

% Random sampling of the state vector
s_meas = [];
time_meas = [];

for i = 1:size(S,1)
    s_meas = [s_meas; S{i}(:,2:end)];
    time_meas = [time_meas; S{i}(:,1)];
end

s_meas = reshape(s_meas.',[size(s_meas,2)*size(s_meas,1) 1]); 
[time_meas, index] = sort(time_meas); 
s_meas = s_meas(index,:);

PD = 0.98;
index = logical(randsrc(length(time_meas),1,[0, 1; 1-PD, PD]));
time_meas = time_meas(index);
s_meas = s_meas(index,:);

% Measurements
sigma_m = 1e-2;
x = normrnd(s_meas(1)*(1-s_meas(2)^2)*cos(s_meas(end)),sigma_m,size(s_meas,1),1);
y = normrnd(s_meas(1)*(1-s_meas(2)^2)*sin(s_meas(end)),sigma_m,size(s_meas,1),1);

R = [x y]; 
meas = [time_meas R];

% Clutter generation 
Pc = 0.0;                   % Probability of false measurements
Vc = 10;                    % Number of false measurements per orbit 

index = logical(randsrc(length(t), 1, [0, 1; 1-Pc, Pc]));
theta_c = 2*pi*rand(length(t),1);
x = normrnd(cos(theta_c),sigma_m,size(theta_c,1),1);
y = normrnd(sin(theta_c),sigma_m,size(theta_c,1),1);
clutter = [x y];
clutter = [t(index) clutter(index,:)];

if (size(clutter,1) > Vc)
    index = randi([0 size(clutter,1)], Vc,1); 
    clutter = clutter(index,:);
end

% Complete measurement set 
meas = [meas; clutter];
[~, index] = sort(meas(:,1)); 
meas = meas(index,:);

%% Particle definition
% J = 100; 
% Jmax = 100; 
% Mean = linspace(0,2*pi,J).'; 
% Sigma = 1*ones(J,1);
% PHD = PHDFilter(J, Jmax, Mean, Sigma, PS, PD);
% 
% PHD = PHD.DefineDomain(0,2*pi,1000);
% PHD = PHD.birth(1, pi, (0.5*pi)^2);
% 
% PHD = PHD.DefinePruning(1e-8, 4);
% 
% PHD = PHD.AssignLikelihood(@(y,z,P)likelihood_function(y,z,P));

%% Estimation 
% Estimation
% [f, X, N(2,:)] = PF.estimator(tspan, meas, UKF_estimator);

% State estimation 
if (~isempty(X{end}))
    Se = sort(X{end}(1,:));
else
    Se = [];
end
St = [];
for i = 1:size(S,1)
    if (S{i}(end,1) >= tspan(end))
        St = [St; S{i}(end,2)];
    end
end
St = mod(St(:,end),2*pi);
St = sort(St.');

%% Results 
figure
hold on
plot(theta, f(:,end), 'r')
for i = 1:size(S,1)
    if (S{i}(end,1) >= tspan(end))
        xline(mod(S{i}(end,2), 2*pi),'k'); 
    end
end
xlabel('$\theta$')
ylabel('$p(\theta)$')
grid on;

x = r.*cos(theta);
y = r.*sin(theta);

figure
view(3)
hold on
plot3(r.*cos(theta),r.*sin(theta),zeros(1,length(theta)))
plot3(x,y,f(:,end))
stem3(r*cos(St.'), r*sin(St.'), 10*ones(size(St)), 'k');
stem3(r.*cos(Se),r.*sin(Se), 10*ones(size(Se)), '*')
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
    plot3(r.*cos(theta),r.*sin(theta),zeros(1,length(theta)));

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
        V = scatter3(r.*cos(se),r.*sin(se), 10*ones(size(se)), '*');

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
% Orbital dynamics
function [ds] = dynamics(mu, t, s)
    ds = [zeros(length(s)-1,1); sqrt(mu/s(1)^3)];
end

% Kinematic proposal 
function [theta_plus, sigma_plus] = kinematic_proposal(r, sigma_r, time_step, theta, sigma)
    % State update
    theta_plus = theta + time_step*r^(-3/2);

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

function set_graphics()
    %Set graphical properties
   set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
    set(groot, 'defaultAxesFontSize', 11); 
    set(groot, 'defaultAxesGridAlpha', 0.3); 
    set(groot, 'defaultAxesLineWidth', 0.75);
    set(groot, 'defaultAxesXMinorTick', 'on');
    set(groot, 'defaultAxesYMinorTick', 'on');
    set(groot, 'defaultFigureRenderer', 'painters');
    set(groot, 'defaultLegendBox', 'off');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultLegendLocation', 'best');
    set(groot, 'defaultLineLineWidth', 1); 
    set(groot, 'defaultLineMarkerSize', 3);
    set(groot, 'defaultTextInterpreter','latex');
end