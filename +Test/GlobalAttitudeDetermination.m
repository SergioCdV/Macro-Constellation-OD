%% In-orbit validation %% 
% Sergio Cuevas del Valle

%% Global attitude determination
% This file provides the function to provide global attitude determination for a given spacecraft.

clear; 
close all; 
clc;
rng(1); 
format long; 

%% Sensor information 
% Covariance matrices
Rm = 0.6^2 * eye(3);                % Magnetometer covariance matrix
Rs = 1e-2^2 * eye(3);               % Sun vector covariance matrix
Rg = deg2rad(0.07)^2 * eye(3);      % Gyroscope covariance matrix

%% Telemetry reading
% Epoch
t = (0:5:3500).'; 
D = datetime([2022 06 11]) + seconds(t);
JD = juliandate(D);
d = decyear(D);

UTC = datetime(JD, 'Format','yyyy MM dd HH mm ss.SSS', 'ConvertFrom', 'juliandate');
UTC.TimeZone = 'Z';
UTC = [year(UTC) month(UTC) day(UTC) hour(UTC) minute(UTC) second(UTC)];

% Spacecraft COE
s0 = [6900e3 1e-3 0 deg2rad(97) 0 0];

%% Pre-processing of the measurements 
% Spacecraft orbit
Earth_obj = Bodies.Earth();
SC_orbit = Orbit(Earth_obj.mu, 'COE', s0, JD(1)).SetFinalEpoch(JD(end)).SetCurrentEpoch(JD(end)).DefineJ2Problem(Earth_obj.J2, Earth_obj.Re).Normalize(true, Earth_obj.Re).AddPropagator('Osculating J2', 1);
SC_orbit = SC_orbit.Propagate().Normalize(false, Earth_obj.Re);

options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);
[~, ASV] = ode113(@(t,s)DiffEvolution(t,s), SC_orbit.StateEvolution(:,1), [0;0;0;1;0;0;0.1], options);

index = zeros(1,length(SC_orbit.StateEvolution(:,1)));
k = 1;
for i = 1:size(SC_orbit.StateEvolution(:,1),1)

    if (SC_orbit.StateEvolution(i,1) >= t(k))
        index(i) = true; 
        k = k + 1;
    end
end

OrbitalStateVector = SC_orbit.ChangeStateFormat('ECI').StateEvolution(logical(index), 2:end);
ASV = ASV(logical(index),:);

% Reference sun vector 
nref(:,1:3) = Astrodynamics.SunEphemeris(JD.').';
nref = nref - OrbitalStateVector(:,1:3);
nref = nref ./ sqrt(dot(nref, nref, 2));

% Transformation to LLA coordinates 
LLA = eci2lla(OrbitalStateVector(:,1:3), UTC);

% Reference geomagnetic field
[nref(:,4:6), ~, ~, ~, ~, ~, ~, ~, ~] = igrfmagm(LLA(:,end), LLA(:,1), LLA(:,2) , d, "13");
nref(:,4:6) = nref(:,4:6) ./ sqrt(dot(nref(:,4:6), nref(:,4:6), 2));

% Body measurements, Synthetic telemetry 
omega = mvnrnd(ASV(:,5:7), Rg, length(t));
nmeas = zeros(length(t), 6); 
for i = 1:length(t)
    y1 = QuaternionAlgebra.RotateVector(ASV(i,1:4).', nref(i,1:3).');
    y2 = QuaternionAlgebra.RotateVector(ASV(i,1:4).', nref(i,4:6).');
    nmeas(i,:) = [y1(1:3,1).' ./ norm(y1(1:3,1).'), y2(1:3,1).' ./ norm(y2(1:3,1).')];
end

%% Processing of the measurements 
Measurements(:,1) = num2cell(t);

for i = 1:length(t)
    % Determine the type of measurement 
    if ( all( ~isnan([nmeas(i,:) omega(i,:)]) ) )
        measurement_type = 'ALL'; 
    elseif ( any( isnan(nmeas(i,1:3)) ) && any( isnan(nmeas(i,4:6)) ) )
        measurement_type = 'GYRO';
    elseif ( any( isnan(nmeas(i,4:6)) ) )
        measurement_type = 'SUN';
    else
        measurement_type = 'MAG';
    end

    % Generate the measurement cell
    switch (measurement_type)
        case 'SUN'
        % Compute the reference vector
        r = nref(i,1:3).';
        r = r ./ norm(r);

        % Normalize the measurement vector
        y = nmeas(i,1:3).' ./ norm(nmeas(i,1:3).'); 

        % Assemble the cell
        Measurements{i,2} = [y; omega(i,:).'];
        Measurements{i,3} = @(s)MeasurementModel(s, r);
        Measurements{i,4} = blkdiag(Rs, Rg);

        case 'MAGNETIC'  
        % Compute the reference vector
        r = nref(i,4:6).';
        r = r ./ norm(r);

        % Normalize the measurement vector
        y = nmeas(i,4:6).' ./ norm(nmeas(i,4:6).'); 

        % Assemble the cell
        Measurements{i,2} = [y; omega(i,:).'];
        Measurements{i,3} = @(s)MeasurementModel(s, r);
        Measurements{i,4} = blkdiag(Rm, Rg);

        case 'GYRO'
        % Assemble the cell
        Measurements{i,2} = omega(i,:).';
        Measurements{i,3} = @(s)GyroModel(s);
        Measurements{i,4} = Rg;

        case 'ALL'
        % Compute the reference vector
        r1 = nref(i,1:3).';
        r1 = r1 ./ norm(r1);

        r2 = nref(i,4:6).';
        r2 = r2 ./ norm(r2);

        r = [r1; r2];

        % Normalize the measurement vector
        y1 = nmeas(i,1:3).' ./ norm(nmeas(i,1:3).');
        y2 = nmeas(i,4:6).' ./ norm(nmeas(i,4:6).');
        y = [y1; y2];

        % Assemble the cell
        Measurements{i,2} = [y; omega(i,:).'];
        Measurements{i,3} = @(s)FullMeasurementModel(s, r);
        Measurements{i,4} = blkdiag(Rs, Rm, Rg);
    end
end

%% Attitude UKF estimation 
% Some constants 
Q = 1e-1 * eye(6);      % Process noise covariance

% Resolution of Wahba's problem as initial conditions for the UKF 
% WSolver = Filters.WahbaSolver();
% [q0, Sigmaq] = WSolver.Davenports(ones(1, 2 * size(nmeas(1,:),1)), [nmeas(1,1:3).' nmeas(1,4:6).'], [nref(1,1:3).' nref(1,4:6).']);

omega0 = omega(1,:).';                      % Initial angular velocity

% Initial conditions 
s0 = [rand(4,1); omega0];                          % Complete state vector
s0(1:4,1) = s0(1:4,1) / norm(s0(1:4,1));
Sigma0 = 1e-6 * eye(6);                            % Initial covariance matrix
Sigma0(1:3,1:3) = 1e-6 * eye(3);                   % Initial MPR covariance

% Propagation model 
SS = @(s, time_step)StateEvolution(s, time_step);

%% Recursive estimation 
% Construction of the estimator 
tspan = t;
StateDim = 6;
USQUE_FC = Filters.USQUE('UKF-A', 2, 1e-3, 0, 1);
USQUE_FC = USQUE_FC.AssignStateProcess(StateDim, SS);
USQUE_FC = USQUE_FC.Init();
USQUE_FC = USQUE_FC.InitConditions(s0, Sigma0);
USQUE_FC.InitFlag = true;

% Preallocation     
S = zeros(7,length(tspan));                            % State estimation
P = zeros(size(Q,1) * size(Q,2), length(tspan));        % Covariance matrix
time = zeros(1,length(tspan));                          % Iteration time
E = zeros(1,length(tspan));                             % Differential entropy

% Recursive Bayes estimation 
last_epoch = tspan(1);                                  % Initial epoch
meas_index = 1;                                         % Measurement index

fprintf('Initialization... \n');
fprintf('------------------------------\n');
fprintf('Running the filter: \n');

% Main loop
i  = 1;
while (i <= length(tspan))
    tic 

    % Check for new measurements to process
    new_measurements = 0;

    if ( ~isempty(Measurements) )
        GoOn = true;
        while (GoOn)
            if (meas_index + new_measurements > size(Measurements,1))
                new_measurements = size(Measurements,1) - meas_index;
                GoOn = false;
            elseif (Measurements{meas_index + new_measurements,1} <= tspan(i))
                new_measurements = new_measurements + 1; 
            else
                GoOn = false;
            end
        end

        % Enable the corrector step
        if ( new_measurements > 0 ) 
            measurement_flag = true;
        else
            measurement_flag = false;
        end
    else
        measurement_flag = false;
    end

    % Perform the correction step if new measurements are available 
    if (measurement_flag)
        for j = 0:new_measurements-1
            % Prepare the measurement process
            prop_epoch = Measurements{meas_index + j, 1};
            z = Measurements{meas_index + j, 2};
            MM = Measurements{meas_index + j, 3};
            R = Measurements{meas_index + j, 4};
            USQUE_FC = USQUE_FC.AssignObservationProcess(size(z,1), MM).AdditiveCovariances(Q, R);

            % Propagation step and weight proposal using the kinematic prior
            [sigma, State, Sigma] = USQUE_FC.PropagationStep(prop_epoch-last_epoch);

            % Correction step 
            [State_C, Sigma_C, ~, ~] = USQUE_FC.CorrectionStep(sigma, State, Sigma, z);

            % New initial conditions 
            State_C = State_C([7:size(State_C,1) 4:6],1);
            USQUE_FC.InitFlag = true;
            USQUE_FC = USQUE_FC.InitConditions(State_C, Sigma_C);
            last_epoch = prop_epoch;
        end

        % Sanity check on the number of processed measurements 
        meas_index = meas_index + new_measurements;
        
    else
        % Propagate to the new epoch the clustered states
        if (last_epoch ~= prop_epoch)
            % Particle propagation
            [~, State_C, Sigma_C] = USQUE_FC.PropagationStep(prop_epoch-last_epoch);
            State_C = State_C([7:size(State_C,1) 4:6],1);
            USQUE_FC = USQUE_FC.InitConditions(State_C, Sigma_C);
        end

        last_epoch = prop_epoch;
    end

    future = i + max(1,new_measurements);

    % Save results
    S(:,i) = State_C;                                   % Predicted state
    P(:,i) = reshape(Sigma_C, [], 1);                   % Predicted covariance
    E(i) = 0.5 * log(det(2*pi * exp(1) * Sigma_C));     % Entropy characterization 

    for j = 1:new_measurements-1
        S(:,i+j) = S(:,i);
        E(i+j) = E(i);
    end

    time(i) = toc;
    fprintf('Iteration running time: %.4f s\n', time(i));

    % Update the timer 
    i = future;
end

fprintf('------------------------------\n');
fprintf('Bayes filter recursion finished.\n');
fprintf('Total running time: %.4f s\n', sum(time));

%% Results
figure
hold on 
plot(tspan, ASV(:,1:4))
scatter(tspan, S(1:4,:))
hold off
xlabel('Mission time [s]')
ylabel('$\hat{\textbf{q}}$')
grid on;

e = zeros(1,size(S,2)); 
for i = 1:size(S,2)
    error = QuaternionAlgebra.right_isoclinic(ASV(i,1:4).') * QuaternionAlgebra.quaternion_inverse(S(1:4,i));
    e(i) = 2 * acos(error(4));%QuaternionAlgebra.log_map(error, [0;0;0;1]);
end

figure
hold on 
plot(tspan, e)
hold off
xlabel('Mission time [s]')
ylabel('$\textbf{e}_q$')
grid on;

figure
hold on 
plot(tspan, rad2deg( ASV(:,5:7)-S(5:7,:).' ) )
hold off
xlabel('Mission time [s]')
ylabel('$\hat{\omega}$ [deg/s]')
grid on;

figure
hold on 
plot(tspan, E)
hold off
xlabel('Mission time [s]')
ylabel('Diff. entropy $E$')
grid on;

%% Auxiliary functions 
% Propagation model for the state 
function [sp] = DiffEvolution(time_step, s)
    % State variables 
    q = s(1:4,:);           % Quaternion
    omega = s(5:7,:);       % Angular velocity

    % Kinematics propagation 
    dq = 0.5 * QuaternionAlgebra.right_isoclinic( [omega; 0] ) * q + 1e-3 * (1-q.'*q) * q;
    sp = [dq; zeros(3,1)];
end

function [sp] = StateEvolution(s, time_step)
    % Constant 
    One = [0;0;0;1];        % Identity quaternion

    % State variables 
    q = s(1:4,:);           % Quaternion
    omega = s(5:7,:);       % Angular velocity

    Omega = omega * time_step/2;

    % Kinematics propagation 
    qp = s(1:4,:); 
    for i = 1:size(s,2)
        qp(:,i) = QuaternionAlgebra.right_isoclinic( QuaternionAlgebra.exp_map(Omega(:,i), One) ) * q(:,i);
    end

    sp = [qp; omega];
end

% Measurement models 
function [y] = MeasurementModel(s, r)
    % Normalize the reference vector 
    r = r/norm(r); 

    % Compute the rotation
    y = zeros(size(r,1), size(s,2));

    for i = 1:size(s,2)
        y(:,i) = QuaternionAlgebra.RotateVector(s(1:4,i), r);
    end

    % Final measurement vector
    y = [y(1:3,:); s(5:7,:)];
end

function [y] = GyroModel(s)
    % Final measurement vector
    y = s(5:7,:);
end

function [y] = FullMeasurementModel(s, r)
    % Compute the rotation
    y = zeros(size(r,1), size(s,2));

    for i = 1:size(s,2)
        y(1:3,i) = QuaternionAlgebra.RotateVector(s(1:4,i), r(1:3,1));
        y(4:6,i) = QuaternionAlgebra.RotateVector(s(1:4,i), r(4:6,1));
    end

    % Final measurement vector
    y = [y; s(5:7,:)];
end