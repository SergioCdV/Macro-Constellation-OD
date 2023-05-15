%% Combinatorial constellation macro-determination 
% Date: 02/02/2023
% Author: Sergio Cuevas del Valle

%% IOD TEST I %%
% This script provides the test for the IOD filter %

close all 
clear 

%% Test input 
tspan = linspace(0,1,1e2);
Measurements = num2cell(tspan).';

for i = 1:size(Measurements,1)
    Measurements{i,2} = @(x)(dot(x,x));
end

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


%% Create the filter 
IOD_filter = Filters.IOD_filter(10, 10, 5, .98, 1);

%% Run the filter
tic
[f, x, n] = IOD_filter.BayesRecursion(tspan, Measurements);
running_time = toc;