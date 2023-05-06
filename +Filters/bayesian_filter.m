%% Combinatorial constellation macro-determination 
% Date: 27/10/2022
% Author: Sergio Cuevas del Valle

%% Toy model for a multi-state Bayesian filter implementation
% The dynamics of a single spacecraft in a equatorial, planar circular
% orbit are analysed, based on inertial measurements of the radial vector

set_graphics();

%% Input data 
mu = 1;                 % Gravitational parameter
r = 10;                 % Radial vector 
theta0 = [0 0];      % Initial conditions
omega = sqrt(mu/r^3);   % Mean motion
T = 1/omega;            % Characteristic time

%% Integration and measurements generation 
% True orbital dynamics
options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
tspan = 0:5e-1:2*pi*T; 

[t, s] = ode45(@(t,s)dynamics(mu,omega,t,s), tspan, theta0, options);

% Random sampling of the state vector
time_index = randi([0 1],length(tspan),2);
time_index = logical(time_index);

s_meas = [s(time_index(:,1),1); s(time_index(:,2),2)]; 
time_meas =  [t(time_index(:,1)); t(time_index(:,2))]; 

[time_meas, re_order] = sort(time_meas); 
s_meas = s_meas(re_order,1);

% Measurements
x = normrnd(r * cos(s_meas),1,size(s_meas,1),1);
y = normrnd(r * sin(s_meas),1,size(s_meas,1),1);
meas = [time_meas x y];

Sigma = (1)^2*eye(2);

%% Estimation 
n = 1e3;                      % Number of particles
f = ones(1,n);                % Uniform distribution                
f(1,:) = f(1,:)/(2*pi);       % Initial angular variable pdf

[theta, f] = estimator(r, omega, meas, Sigma, f);

%% Results 
figure
hold on
stem(theta, f, 'r')
xline(mod(s(end,1),2*pi),'k');
xline(mod(s(end,2),2*pi),'k');
legend('PDF', '$\theta_f$')
xlabel('$\theta$')
ylabel('$p(\theta)$')
grid on;

x = r.*cos(theta);
y = r.*sin(theta);

figure
view(3)
hold on
plot3(r.*cos(theta),r.*sin(theta),zeros(1,length(theta)))
plot3(x,y,f)
stem3(r.*cos(s(end)),r.*sin(s(end)),1,'filled')
zlabel('$p(\theta)$')
grid on;

%% Auxiliary functions 
% Bayesian estimation 
function [thetaf, posterior] = estimator(g, omega, meas, Sigma, f0)
    % Preallocation 
    n = size(f0, 2);
    theta = linspace(0,2*pi,n).'; 

    prior = f0(1,:).';

    % Compute the likelihood function
    [m, J] = measurement_model(g(1), theta);
    pz = likelihood_function(m, meas(1,2:3), Sigma)./J;

    % Compute the Bayesian posterior 
    posterior = pz.*prior;

    % Normalize
    posterior = posterior/sum(posterior);

    % Bayesian estimation
    for i = 2:size(meas,1)
        % Step ahead the prior
        thetaf = step_ahead(theta, omega, meas(i,1));

        % Compute the likelihood function
        [m, J] = measurement_model(g(1), thetaf);
        pz = likelihood_function(m, meas(i,2:3), Sigma)./J;

        % Compute the Bayesian posterior 
        posterior = pz.*prior;

        % Normalize
        alpha = sum(posterior);
        posterior = posterior/alpha;

        % Next step 
        prior = posterior;
    end

    thetaf = mod(thetaf,2*pi);
end

% Orbital dynamics
function [dtheta] = dynamics(mu,omega,t,s)
    dtheta = [omega; 
              omega];
end

% Measurement function 
function [m, J] = measurement_model(rho, theta)
    m = [rho*cos(theta) rho*sin(theta)]; 
    J = rho;
end

% Measurement likelihood function
function [pz] = likelihood_function(z, mu, Sigma) 
    % Compute the Gaussian estimation
    res = z-mu; 

    pz = zeros(size(res,1),1);

    for i = 1:size(res,1)
        pz(i) = 1/sqrt(det(2*pi*Sigma))*exp(-0.5*res(i,:)*Sigma^(-1)*res(i,:).');
    end
end

% STM function 
function [theta] = step_ahead(theta0, omega, dt)
    theta = theta0+dt*omega;
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