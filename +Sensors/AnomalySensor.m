%% Constellation macro-orbit determination 
% Date: 06/02/2023
% Author: Sergio Cuevas del Valle

%% Class implementation of an AnomalySensor object

classdef AnomalySensor< Sensors.AbstractSensor
    properties  
    end 

    methods 
        % Constructor 
        function [obj] = AnomalySensor(myInitialEpoch, myInitialState, mySigma, myPD)
            myStateDim = 1;
            myMeasDim = 1;
            obj = obj@Sensors.AbstractSensor(myStateDim, myMeasDim, myInitialEpoch, myInitialState, mySigma, myPD);
        end

        % Configuration of the clutter model 
        function [obj] = ConfigClutter(myPC, myNC)
            if (myPC >= 0 && myPC <= 1)
                obj.PC = myPC;
            else
                error('No valid detection probability has been input.'); 
            end

            if (myNC > 0)
                obj.NC = myPC; 
            else
                error('The average number of false measurements is not valid.');
            end
        end

        % Observations
        function [timestamp, meas, StateEvolution] = Observe(obj, Orbit, Epoch)
            % Propagate the orbit 
            AuxOrbit = Orbit.SetCurrentEpoch(Epoch).Propagate().ChangeStateFormat('COE');

            % Propagate the observer and take the measurements
            AuxOrbit = AuxOrbit.Normalize(false, 1);
            Tspan = AuxOrbit.StateEvolution(:,1);
            [Tspan, StateEvolution, index] = obj.Dynamics(AuxOrbit.InitialEpoch, obj.State, Tspan); 
            AuxOrbitEvolution = AuxOrbit.StateEvolution(index,2:end);

            if (~isempty(StateEvolution))
                % Restrict the orbit state evolution                    
                [timestamp, meas] = obj.ObservationProcess(Tspan, AuxOrbitEvolution, StateEvolution);
    
                % Apply the probability of detection
                index = logical(randsrc(length(timestamp), 1, [0, 1; 1-obj.PD, obj.PD]));
                timestamp = timestamp(index); 
                meas = meas(index,:); 
                StateEvolution = StateEvolution(index,:);

                % Add some noise
                for i = 1:size(meas,1)
                    meas(i,:) = mvnrnd(meas(i,:), obj.Sigma, 1);
                end

                % Add clutter
                index = logical(randsrc(size(meas,1), 1, [0, 1; 1-obj.PC, obj.PC]));
                clutter = rand(length(index), 3); 
                clutter = clutter./sqrt(dot(clutter,clutter,2));
                clutter = clutter .* ( mean(sqrt(dot(meas(index,:),meas(index,:),2))) );
                clutterTime = timestamp(index,:);
                clutterState = StateEvolution(index,:);
                
                if (size(clutter,1) > round(1.5 * obj.NC))
                    index = randi([0 size(clutter,1)], obj.NC, 1); 
                    clutter = clutter(index,:);
                    clutterTime = clutterTime(index,:);
                    clutterState = clutterState(index,:);
                end

                % Final assembly
                meas = [meas clutter(index,:)];
                timestamp = [timestamp; clutterTime];
                StateEvolution = [StateEvolution; clutterState];
                [timestamp, index] = sort(timestamp);
                meas = meas(index,:);
                StateEvolution = StateEvolution(index,:);

            else
                timestamp = []; 
                meas = []; 
                StateEvolution = [];
            end
        end

        % Observation process 
        function [t, meas] = ObservationProcess(obj, Tspan, Orbit, StateEvolution)
            % Preallocaton
            meas = []; 
            t = []; 

            % Observation
            for i = 1:length(Tspan)         
                meas = [meas; Orbit(i,6)];
                t = [t; Tspan(i)];
            end
        end

        % Likelihood function
        function [q] = LikelihoodFunction(obj, Sigma, z, y)
            res = mod(y,2*pi)-mod(z,2*pi);
            q = exp((-0.5*res.'*Sigma^(-1)*res))/sqrt(det(Sigma)*(2*pi)^(size(Sigma,1)));
        end
    end
    
    methods (Access = private)
        % Dynamics 
        function [epoch, StateEvolution, index] = Dynamics(obj, Epoch, State, Tspan)
            % Check if an initial state has been given 
            if (~isempty(obj.State))
                % Check the tspan 
                index = Epoch + Tspan >= obj.InitialEpoch;
                epoch = Epoch + Tspan(index) / (24 * 3600);
    
                % Copy the state
                StateEvolution = repmat(State, length(Tspan(index)), 1);
                Tspan = Tspan(index);
            else
                error('No obsever state has been given.');
            end
        end
    end
end