%% Constellation macro-determination %%
% Date: 06/12/2023

%% Constellation orbit determination. Scenario IOD %%
% This script provides the IOD scenario in which
% to test the on-manifold PHD to identify individual targets %

close all 
clear 
rng(1);         % Reproducibility

Generation = false;

if (Generation)
    SciTech2024.GenMeasurements;
else
    load SciTech2024_IOD.mat;
end

%% Estimation: IOD
% Estimator configuration
PD = 0.9;                                               % Probability of detection
PHIOD_filter = Filters.PHIOD_filter(1, 1e3, PS, PD);    % Filter initialization

PHIOD_filter.Jmax = 100;                                % Maximum number of mixture components

PHIOD_filter.epsilon = 0;                               % Perturbation parameter

PHIOD_filter.Lmin = 1.03;                               % Box constraint
PHIOD_filter.Lmax = 2;
PHIOD_filter.emax = 0.2;

% Estimation
[X, N_hat, Prior, E_IOD] = PHIOD_filter.BayesRecursion(ObservationSpan, Measurements);

% Plane states
[Planes] = PHIOD_filter.StateEstimation(X{end}(1:end-1,:), X{end}(end,:));
    