%% Macro-orbit determination 
% Keplerian uncertainty propagation 
% Date 08/07/2023

clear; 
clc; 
close all; 
set_graphics();

%% Input data 
% Meshgridding 
a = linspace(1.03, 2, 100);         % Semimajor axis grid [Earth radii]
omega = linspace(-0.04, 0.04, 100);         % Time element [rad]

% Define the initial PDF parametes 
Mu = [1.4322; 0];
Q = diag([0.25^2 0.02^2]); 
invQ = Q^(-1);

% Propagation time 
pepochs = linspace(0, 1, 50); 

% Preallocation for speed 
f = cell(length(pepochs),1);
M = cell(length(pepochs),1);
u = zeros(length(a), length(omega));

for i = 1:length(pepochs)
    M{i} = omega + sqrt(1/a(1)^3) * pepochs(i);
    for j = 1:length(a)
        % Computation of the pdf 
        Mm = M{i} - sqrt(1/a(j)^3) * pepochs(i);
        x = [repmat(a(j), 1, length(omega)); Mm];
        u(j,:) = exp(-0.5 * dot((x-Mu), invQ * (x-Mu))) / (2*pi*det(Q));
    end

    % Save the probabily function
    f{i} = u;
end

% Plot results
filename = ['phase_uncertainty.gif'];
h = figure(1);

ylim([0.2 0.6])
grid on;
for i = 1:length(pepochs)
    % Create the meshgrid 
    [A, o] = meshgrid(a, M{i});

    % Phase space evolution 
    contourf(A, o, f{i});
    view([0 90])
    xlabel('$a$ [$R_e$]')
    ylabel('$M$')
    zlabel('$\rho$')
    colorbar
    pause(2)

    if (true)
        drawnow;
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        if (i == 1) 
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1e-3); 
        else 
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1e-3); 
        end 
    end
end

%% Set up cool graphics 
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
