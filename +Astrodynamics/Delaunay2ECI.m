

function [s] = Delaunay2ECI(D)
    % Delaunay set to the Cartesian state transformation
    M = D(1);           % Mean anomaly
    omega = D(2);       % AoP
    Omega = D(3);       % RAAN
    L = D(4);           % Delaunay action
    G = D(5);           % Angular momentum
    H = D(end);         % Polar angular momentum

    % Keplerian functions
    e = real(sqrt(1-(G/L)^2));
    cos_i = H/G; 

    % Perifocal quaternion
    qp(1,1) = sin(acos(cos_i)/2) * cos((Omega-omega)/2);
    qp(2,1) = sin(acos(cos_i)/2) * sin((Omega-omega)/2);
    qp(3,1) = cos(acos(cos_i)/2) * sin((Omega+omega)/2);
    qp(4,1) = cos(acos(cos_i)/2) * cos((Omega+omega)/2);

    % Solve for the true anomaly
    nu = Astrodynamics.KeplerSolver(e,M);

    % Compute the radial and velocity distance in the perifocal frame
    R = G^2./(1+e*cos(nu));
    r = R .* [cos(nu); sin(nu); zeros(1,length(nu))];
    v = [-sin(nu); e + cos(nu); zeros(1,length(nu))] / G;

    Q = QuaternionAlgebra.right_isoclinic( QuaternionAlgebra.quaternion_inverse(qp) );

    % Compute the ECI osculating coordinates
    aux(:,1) = QuaternionAlgebra.right_isoclinic( [r; 0] ) * qp;
    aux(:,2) = QuaternionAlgebra.right_isoclinic( [v; 0] ) * qp;
    aux = Q * aux; 
    State = [aux(1:3,1); aux(1:3,2)];

    s = real(State);
end