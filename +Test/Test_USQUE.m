

%% USQUE TEST 


clear; 
clc; 

%% Input data 
t = 0:1e-3:10;
omega0 = zeros(length(t),3);
omega = omega0 + 1e-3 * rand(length(t),3);
qt = repmat([0;0;0;1],length(t),1);

r = [0;1;0];
m = r.' + 1e-9 * rand(length(t),3);
m = m./sqrt(dot(m,m,2));

Q = 1e-7 * eye(6);
R = 1e-7 * eye(6);

% Initial conditions 
q0 = [sqrt(2)/2; 0; 0; sqrt(2)/2; ones(3,1)];
Sigma0 = 1e-1 * eye(6);

%% Recursive estimation 
SS = @(s,dT)StateEvolution(s,dT);

% Create the filter 
a = 2;
USQUE_ex = Filters.USQUE('UKF-S', 2, 1e-1, 0, 1).AdditiveCovariances(Q, R);
USQUE_ex = USQUE_ex.AssignStateProcess(6, SS);
USQUE_ex = USQUE_ex.Init();
USQUE_ex = USQUE_ex.InitConditions(q0, Sigma0);

tic
for i = 1:length(t)
    MM = @(s)MeasurementModel(s, r);
    USQUE_ex = USQUE_ex.AssignObservationProcess(6, MM);

    % Propagation 
    [sigma, State, Sigma] = USQUE_ex.PropagationStep(t(i));

    % Correction
    [Statec(:,i), Sigmac, ~, ~] = USQUE_ex.CorrectionStep(sigma, State, Sigma, [m(i,:) omega(i,:)].');

    USQUE_ex = USQUE_ex.InitConditions([Statec(7:end,i); Statec(4:6,i)], Sigmac);

    % Eigenvalues monitoring 
    [~, V] = eig(Sigma*Sigma.'); 
    lambda(i) = min(diag(V));
end
toc

%% Auxiliary functions 
function [qp] = StateEvolution(s, dT)
    qp = s; 
    for i = 1:size(s,2)
        qp(1:4,i) = QuaternionAlgebra.right_isoclinic(s(1:4,i)) * QuaternionAlgebra.exp_map([s(5:7,i)*dT/2; 0], [0;0;0;1]);
    end
end

function [y] = MeasurementModel(q, r)
    for i = 1:size(q,2)
        y(:,i) = QuaternionAlgebra.left_isoclinic(q(1:4,i)).' * QuaternionAlgebra.right_isoclinic(q(:,i)) * [r; 0];
    end
    y = [y(1:3,:); q(5:7,:)];

end