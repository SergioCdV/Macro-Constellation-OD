

function [sigma, State, Sigma] = PropagationStep(obj, time_step)
    % Propagate the noise 
%     Q = time_step / 2 * [(obj.Q(1,1)-time_step^2/6*obj.Q(2,2)) * eye(3) zeros(3); zeros(3) obj.Q(2,2) * eye(3)];
%     obj.Q = Q;

    % Generate sigma points 
    sigma = obj.sigma_points([zeros(3,1); obj.State(5:7,1)], obj.Sigma, obj.Q);

    % Generation of the quaternions 
    X = [QuaternionAlgebra.MPR2Quat(obj.a, obj.f, sigma(1:3,:), true); sigma(4:6,:)];

    for i = 1:size(X,2)
        X(1:4,i) = QuaternionAlgebra.right_isoclinic(X(1:4,i)) * obj.State(1:4,1);
    end

    % Propagation of sigma points 
    X = feval(obj.StateModel, X, time_step);

    % Recover the MPR 
    dq = X(1:4,:); 
    Q = QuaternionAlgebra.quaternion_inverse(X(1:4,1));
    for i = 1:size(dq,2)
        dq(:,i) = QuaternionAlgebra.right_isoclinic(X(1:4,i)) * Q;
    end

    sigma = [QuaternionAlgebra.MPR2Quat(obj.a, obj.f, dq(1:4,:), false); X(5:7,:)];

    % State and covariance prediction 
    switch (obj.Algorithm)
        case 'UKF-A'
            [State, Sigma] = UKFA_prediction(obj, sigma);
        case 'UKF-S'
            [State, Sigma] = UKFS_prediction(obj, sigma);
    end

    % Quaternion parts
    dq = QuaternionAlgebra.MPR2Quat(obj.a, obj.f, State(1:3,:), true);
    QuatState = QuaternionAlgebra.right_isoclinic(dq) * obj.State(1:4,1);
    State = [State; QuatState];
    sigma = [sigma; X(1:4,:)];
end