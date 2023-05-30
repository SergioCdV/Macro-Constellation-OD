
% Transformation from Delaunay to Cartesian elements
function [State] = ParticleState(obj, SensorModality, particle, nu)
    for j = 1:length(nu)
        % Delaunay elements of the particle 
        qp = particle(1:4,1);
        L = particle(5,1);
        G = particle(6,1);
        H = particle(7,1);
        M = nu(j);
     
        % Compute the RAAN and AoP from qp 
        diff = atan2(qp(2,1), qp(1,1));
        plus = atan2(qp(3,1), qp(4,1));
        Omega = 2 * (plus+diff);
        omega = 2 * plus - Omega;
    
        % Assemble the set
        D = [M omega Omega L G H].';   
    
        switch (SensorModality)
            case 'ANOMALY'
                State(:,j) = [zeros(1,5) M].';
    
            case 'DELAUNAY'
                state = Astrodynamics.Delaunay2COE(1, D, true);
                state(1) = state(1) * obj.Re;
                State(:,j) = state.';
    
            otherwise
%                 Do = Astrodynamics.Brouwer_solution(obj.epsilon, D);
%                 State(:,i) = Astrodynamics.Delaunay2ECI(Do);
%                 State(1:3,i) = obj.Re * State(1:3,i);
%                 State(4:6,i) = obj.Re/obj.Tc * State(4:6,i);

                State(:,j) = Astrodynamics.Lara_solution(obj.epsilon, D);
                State(1:3,j) = obj.Re * State(1:3,j);
                State(4:6,j) = obj.Re/obj.Tc * State(4:6,j);
        end    
    end
end