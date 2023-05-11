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


%% Create the filter 
IOD_filter = Filters.IOD_filter(10, 10, 5, .98, 1);

%% Run the filter
tic
[f, x, n] = IOD_filter.BayesRecursion(tspan, Measurements);
running_time = toc;