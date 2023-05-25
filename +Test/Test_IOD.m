%% Combinatorial constellation macro-determination 
% Date: 02/02/2023
% Author: Sergio Cuevas del Valle

%% IOD TEST I %%
% This script provides the test for the IOD filter %

close all 
clear 

%% Test input 
r0 = 6900e3;                % Characteristic distance of the Earth orbit  
mu = 3.986e14;              % Gravitional parameter of the Earth
Re = 6378e3;                % Reference Earth radius
Tc = sqrt(Re^2/mu);         % Characteristic time
J2 = 1.08263e-3;            % Earth's J2 parameter

%% Some orbital dynamics tests 
X = [1;0;0]; 

Omega = deg2rad( 20 ); 
omega = deg2rad( 0 ); 
i = deg2rad( 90 );

D(1,1) = sin(i/2) * cos((Omega-omega)/2);
D(2,1) = sin(i/2) * sin((Omega-omega)/2);
D(3,1) = cos(i/2) * sin((Omega+omega)/2);
D(4,1) = cos(i/2) * cos((Omega+omega)/2);

x = QuaternionAlgebra.right_isoclinic(D) * (QuaternionAlgebra.right_isoclinic([X; 0]) * QuaternionAlgebra.quaternion_inverse(D));

R1 = [cos(Omega) sin(Omega) 0; -sin(Omega) cos(Omega) 0; 0 0 1];
R2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R3 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];

R = R3*R2*R1;

x2 = R*X;

%% Lara solution test
InitialEpoch = juliandate(datetime('now'));         % Initial epoch in JD
T = 0.5;                                            % Number of days 
EndEpoch = juliandate(datetime('now')+days(T));     % End epoch
Step = 600;                                         % Integration step in seconds
tspan = 0:Step:T * 86400;                           % Relative lifetime in seconds

ElementType = 'COE'; 
ElementSet = [r0 1e-3 deg2rad(0) deg2rad(45) deg2rad(0) 2*pi*rand()];

% Add the orbit to the constellation
AuxOrbit = Orbit(mu, ElementType, ElementSet, InitialEpoch);
AuxOrbit = AuxOrbit.SetCurrentEpoch(EndEpoch); 
AuxOrbit = AuxOrbit.SetFinalEpoch(EndEpoch);
AuxOrbit = AuxOrbit.Normalize(true, Re); 
AuxOrbit = AuxOrbit.ChangeStateFormat('COE');
AuxOrbit = AuxOrbit.DefineJ2Problem(J2, Re);
AuxOrbit = AuxOrbit.AddPropagator('Osculating J2', Step/86400);

AuxOrbit = AuxOrbit.Propagate();
hold on;
AuxOrbit.PlotTrajectory(figure(1), AuxOrbit.InitialEpoch, AuxOrbit.PropagatedEpoch);
%%
L = sqrt(ElementSet(1)/Re);
G = L * sqrt(1-ElementSet(2)^2);
H = G * cos(ElementSet(4));
D = [ElementSet(end) ElementSet(5) ElementSet(3) L G H] + [tspan.'/L^3 zeros(length(tspan),5)];

for i = 1:length(tspan)
    s(i,:) = Astrodynamics.Lara_solution(-J2,D(i,:));
end
scatter3(s(:,1),s(:,2),s(:,3));

%% Filter test
tspan = linspace(0,1,1e2);
Measurements = num2cell(tspan).';

for i = 1:size(Measurements,1)
    Measurements{i,2} = @(x)(dot(x,x));
end

IOD_filter = Filters.IOD_filter(10, 10, 5, .98, 1);

%% Run the filter
tic
[f, x, n] = IOD_filter.BayesRecursion(tspan, Measurements);
running_time = toc;