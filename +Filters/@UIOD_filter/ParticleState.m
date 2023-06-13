
% Transformation from Delaunay to Cartesian elements
function [State] = ParticleState(obj, SensorModality, particle, nu)
    for i = 1:size(particle,2)
        % Delaunay elements of the particle 
        D = Astrodynamics.Delaunay2MyElements(particle(:,i), false);
        D(1,1) = nu;
    
        switch (SensorModality)
            case 'ANOMALY'
                State(:,i) = [zeros(1,5) M].';
    
            case 'DELAUNAY'
                state = Astrodynamics.Delaunay2COE(1, D, true);
                state(1) = state(1) * obj.Re;
                State(:,i) = state.';
    
            otherwise
                State(:,i) = Astrodynamics.Lara_solution(obj.epsilon, D);
                State(1:3,i) = obj.Re * State(1:3,i);
                State(4:6,i) = obj.Re/obj.Tc * State(4:6,i);
        end    
    end
end