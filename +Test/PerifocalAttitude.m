% Macro-orbit determination 
% Keplerian uncertainty propagation 
% Date 08/07/2023

clear; 
clc; 
close all; 
set_graphics();

%% Input data 
J2 = -1.08e-3;

% Polar Delaunay set 
D = [1.05 1e-3 0 pi/4 0 0];
eta = sqrt(1-D(2)^2);
D = Astrodynamics.Delaunay2COE(1, D, false);
tspan = linspace(0, 31*24*2*pi, 1e2); 

omega = zeros(1,6);
omega(3) = (3*J2)/(2*D(4)^7*eta^4)*D(6)/D(5);
omega(2) = (3*J2)/(4*D(4)^7*eta^4)*(1-5*D(6)^2/D(5)^2);
omega(1) = 1/D(4)^3 + (3*J2)/(4*D(4)^7*eta^4)*(1-3*D(6)^2/D(5)^2);
for i = 1:length(tspan)
    D(i,:) = D(1,:) + omega * tspan(i);
end

% Compute the quaternion
C = zeros(length(tspan), 8);
for i = 1:length(tspan)
    qp = [0;0;sin(D(i,1)/2);cos(D(i,1)/2)];
    C(i,:) = Astrodynamics.Delaunay2MyElements(D(i,:).', true).';
    C(i,1:4) = (QuaternionAlgebra.right_isoclinic(qp) * C(i,1:4).').';
end

%% Plots
filename = 'perifocal_evolution.gif';
h = figure(1);
view(3)
hold on
no = [eye(3); zeros(1,3)];
quiver3(zeros(3,1), zeros(3,1), zeros(3,1), no(1:3,1), no(1:3,2), no(1:3,3), 'b', 'LineWidth', 1.0); 

m = 100;
[aa, bb, cc] = sphere(m);
g = surf(aa, bb, cc);
set(g, 'FaceAlpha', 0.1)
shading interp
grid on;
xlabel('$X$')
ylabel('$Y$')
zlabel('$Z$')
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

for i = 1:size(C,1)
    aux = QuaternionAlgebra.right_isoclinic([1 0 0 0].') * QuaternionAlgebra.quaternion_inverse(C(i,1:4).');
    S = QuaternionAlgebra.right_isoclinic(C(i,1:4).') * aux;
    %quiver3(0, 0, 0, S(1,1), S(2,1), S(3,1), 'r', 'LineWidth', 1.0); 
    scatter3(S(1,1), S(2,1), S(3,1), 'r', 'filled'); 
    pause(1)

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
