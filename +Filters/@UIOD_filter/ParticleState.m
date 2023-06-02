
% Transformation from Delaunay to Cartesian elements
function [State] = ParticleState(obj, SensorModality, particle, nu)
    for i = 1:size(particle,2)
        % Delaunay elements of the particle 
        qp = particle(1:4,i);
        L = particle(5,i);
        G = particle(6,i);
        H = particle(7,i);
    
        % Compute the RAAN and AoP from qp 
        diff = atan2(qp(2,1), qp(1,1));
        plus = atan2(qp(3,1), qp(4,1));
        Omega = 2 * (plus+diff);
        omega = 2 * plus - Omega;
        M = nu;
    
        % Assemble the set
        D = [M omega Omega L G H].';   
    
        switch (SensorModality)
            case 'ANOMALY'
                State(:,i) = [zeros(1,5) M].';
    
            case 'DELAUNAY'
                state = Astrodynamics.Delaunay2COE(1, D, true);
                state(1) = state(1) * obj.Re;
                State(:,i) = state.';
    
            otherwise
%                 Do = Astrodynamics.Brouwer_solution(obj.epsilon, D);
%                 State(:,i) = Astrodynamics.Delaunay2ECI(Do);
%                 State(1:3,i) = obj.Re * State(1:3,i);
%                 State(4:6,i) = obj.Re/obj.Tc * State(4:6,i);

                State(:,i) = Astrodynamics.Lara_solution(obj.epsilon, D);
                State(1:3,i) = obj.Re * State(1:3,i);
                State(4:6,i) = obj.Re/obj.Tc * State(4:6,i);
        end    
    end
end