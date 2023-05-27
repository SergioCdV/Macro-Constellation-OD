
% Transformation from Delaunay to Cartesian elements
function [State] = ParticleState(obj, particle)
    % Delaunay elements of the particle 
    qp = particle(1:4,1);
    L = particle(5,1);
    G = particle(6,1);
    H = particle(7,1);
    M = particle(8,1);
 
    % Compute the RAAN and AoP from qp 
    diff = atan2(qp(2,1), qp(1,1));
    plus = atan2(qp(3,1), qp(4,1));
    Omega = 2 * (plus+diff);
    omega = 2 * plus - Omega;

    % Assemble the set
    D = [M omega Omega L G H].';   
    State = [zeros(1,5) D(1)].';

    % Preallocation 
%     Do = Astrodynamics.Brouwer_solution(obj.epsilon, D);
%     State = zeros(6,size(D,2)); 
%     for i = 1:size(D,2)
%         State(:,i) = Astrodynamics.Delaunay2ECI(Do);
%         State(1:3,i) = obj.Re * State(1:3,i);
%         State(4:6,i) = obj.Re/obj.Tc * State(4:6,i);
%     end

%     State = zeros(6,size(D,2)); 
%     for i = 1:size(D,2)
%         State(:,i) = Astrodynamics.Lara_solution(obj.epsilon, D);
%         State(1:3,i) = obj.Re * State(1:3,i);
%         State(4:6,i) = obj.Re/obj.Tc * State(4:6,i);
%     end
end