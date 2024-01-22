%% Plane Identification 
% Date: 13/12/2023

% This script provides several test cases for the proposed membership
% function between spacecraft states (through Milankovitch elements) and
% orbital planes 

%% Input data 
% Generate a pair of measurements 
mu = 1; 

% First particle
e1 = [0; 0; 0.01]; 
h1 = [1; 0; 0];
P1 = 1e-6 * eye(6);

% Second particle 
e2 = [0; 0; 0.01];
h2 = [1; 0; 0];

P2 = 1E-6 * eye(6);
P2(1,2) = 1E-6;
% P2(2,3) = 1E-4;
% P2 = 0.5 * (P2 + P2.');

% Check if the two particles lie in the same orbit 
p = checkOrbit(mu, h1, e1, P1, h2, e2, P2);

%% Auxiliary function 
% Function to check if two states belong to the same plan
function [p] = checkPlane(h1, e1, P1, h2, e2, P2)
    % Tranform to the two S frames 
    [V1, P1] = eig(P1);     % Eigenspectrum of the first covariance
    [V2, P2] = eig(P2);     % Eigenspectrum of the second covariance 

    h1 = V1 * [h1; e1]; 
    h2 = V2 * [h2; e2];

    % Uncorrelated vector measurements
    e1 = h1(4:6,1);
    h1 = h1(1:3,1); 

    e2 = h2(4:6,1);
    h2 = h2(1:3,1); 

    % Compute the error vectors 
    e_h = h2 - h1;      % Error in the angular momentum 
    e_e = e2 - e1;      % Error in the eccentricity vector 
end

% Function to check if two states belong to the same orbit
function [p] = checkOrbit(mu, h1, e1, P1, h2, e2, P2)
    % Tranform to the two S frames 
    [V1, P1] = eig(P1);     % Eigenspectrum of the first covariance
    [V2, P2] = eig(P2);     % Eigenspectrum of the second covariance 

    meas(:,1) = V1 * [h1; e1]; 
    meas(:,2) = V2 * [h2; e2];

    % Mean observations
    h_norm = sqrt(dot( meas(1:3,:), meas(1:3,:), 1));   % Estimated mean angular momentum
    e_norm = sqrt(dot( meas(4:6,:), meas(4:6,:), 1));   % Estimated mean eccentricities

    m = - 0.5 * (1 - e_norm.^2) * mu^2 ./ h_norm.^2;    % Mean energy observations for the orbit (non-singular except for rectilinear trajectories)

    if (sign(m(1)) == sign(m(2)))
        % Transform the covariances 
        if (e_norm(1))
            dA = [meas(1:3,1).'/norm(meas(1:3,1)) zeros(1,3); zeros(1,3) meas(4:6,1).'/norm(meas(4:6,1))];
        else
            dA = [meas(1:3,1).'/norm(meas(1:3,1)) zeros(1,3); zeros(1,6)];
        end
    
        dA = -0.5 * [-2*(1 - e_norm(1)^2)*mu^2/h_norm(1)^3, -2*e_norm(1)*mu^2/h_norm(1)^2] * dA;
        P(1) = dA * P1 * dA.';
    
        if (e_norm(2))
            dA = [meas(1:3,2).'/norm(meas(1:3,2)) zeros(1,3); zeros(1,3) meas(4:6,2).'/norm(meas(4:6,2))];
        else
            dA = [meas(1:3,2).'/norm(meas(1:3,2)) zeros(1,3); zeros(1,6)];
        end
    
        dA = -0.5 * [-2*(1 - e_norm(2)^2)*mu^2/h_norm(2)^3, -2*e_norm(2)*mu^2/h_norm(2)^2] * dA;
        P(2) = dA * P2 * dA.';

        [P, index] = sort(P);
        m = m(index);
    
        % Compute the intersection region between the two Gaussians 
        A = 1 / P(2) - 1 / P(1);
        B = 2 * (-m(2)/P(2) + m(1)/P(1));
        Delta = 4 / (P(1)*P(2)) * ((m(1)-m(2))^2 + (P(1)-P(2)) * log(P(1)/P(2)));
    
        if (A ~= 0)
            X(2) = (-B - sqrt(Delta)) / (2 * A);
            X(1) = (-B + sqrt(Delta)) / (2 * A);
            X = sort(X);
    
            % Compute the probability of both random variables lying in the intersection area 
            P = sqrt(P);
            p(1) = +0.5 * erf( (X(2) - m(1)) / sqrt(P(1) * 2) ) - 0.5 * erf( (X(1) - m(1)) / sqrt(P(1) * 2) );
            p(2) = +0.5 * erf( (X(2) - m(2)) / sqrt(P(2) * 2) ) - 0.5 * erf( (X(1) - m(2)) / sqrt(P(2) * 2) );

            p = max(0, sum(p) * 0.5);
        else
            p = 1;
        end
    else
        p = 0;
    end
end